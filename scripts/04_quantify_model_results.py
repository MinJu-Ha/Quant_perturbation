import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import yaml
from sklearn.linear_model import LogisticRegression
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler


def load_config(path):
    with open(path, "r") as f:
        return yaml.safe_load(f)


def read_expression_matrix(path):
    matrix = pd.read_csv(path, index_col=0)
    matrix.index = matrix.index.astype(str)
    matrix.columns = matrix.columns.astype(str)
    return matrix


def read_chip_targets(path):
    chip = pd.read_csv(path)
    first_col = chip.columns[0]
    return pd.Index(chip[first_col].dropna().astype(str).str.replace("\ufeff", "", regex=False).unique())


def load_metadata(config, expr):
    meta = pd.read_csv(config["paths"]["exported_meta_csv"], index_col=0)
    if "cell_id" in meta.columns:
        meta.index = meta["cell_id"].astype(str)
    meta.index = meta.index.astype(str)
    common = meta.index.intersection(expr.columns)
    if len(common) == 0:
        raise ValueError("No overlapping cell IDs between metadata and expression matrix.")
    return meta.loc[common].copy(), expr.loc[:, common].copy()


def selected_genes(config, expr):
    out = Path(config["paths"]["tables_dir"]) / "selected_bcl6_target_genes.csv"
    if out.exists():
        genes = pd.read_csv(out).iloc[:, 0].dropna().astype(str)
        genes = [g for g in genes if g in expr.index]
        if genes:
            return genes

    chip_targets = read_chip_targets(config["paths"]["chip_targets_csv"])
    expr_upper = {g.upper(): g for g in expr.index}
    targets = []
    for gene in chip_targets:
        hit = expr_upper.get(gene.upper())
        if hit:
            targets.append(hit)
    ko = str(config["perturbation"].get("gene", config["perturbation"]["target_gene_candidates"][0])).upper()
    targets = pd.Index(pd.unique([g for g in targets if g.upper() != ko]))
    if len(targets) < int(config["features"]["min_target_genes"]):
        raise ValueError(f"Too few BCL6 target genes found: {len(targets)}")

    max_targets = int(config["features"]["max_target_genes"])
    if len(targets) > max_targets:
        targets = expr.loc[targets].var(axis=1).sort_values(ascending=False).head(max_targets).index
    pd.Series(targets, name="selected_bcl6_target_gene").to_csv(out, index=False)
    return list(targets)


def fit_classifier(config, meta, expr, genes):
    condition_col = config["metadata"]["condition_column"]
    recovered = set(config["metadata"]["recovered_labels"])
    nonrecovered = set(config["metadata"]["nonrecovered_labels"])
    keep = meta[condition_col].isin(recovered.union(nonrecovered))

    X = expr.loc[genes, meta.index[keep]].T
    y = meta.loc[keep, condition_col].isin(recovered).astype(int).values
    if len(np.unique(y)) < 2:
        raise ValueError("Need both recovered and non-recovered labels for classifier.")

    pipeline = Pipeline([
        ("scaler", StandardScaler(with_mean=True, with_std=True)),
        ("clf", LogisticRegression(
            max_iter=int(config["classifier"]["max_iter"]),
            class_weight="balanced",
        )),
    ])
    pipeline.fit(X, y)
    return pipeline


def distance_to_centroid(matrix_cells_by_genes, centroid):
    return np.linalg.norm(matrix_cells_by_genes - centroid[None, :], axis=1)


def build_perturbed_matrix(model_output, original_expr, genes, cells, output_format):
    output = model_output.reindex(index=genes, columns=cells)
    original = original_expr.loc[genes, cells]
    if output_format == "delta_expression":
        output = output.fillna(0.0)
        return original + output
    output = output.combine_first(original)
    return output


def score_model(config, model, meta, original_expr, perturbed_expr, genes, pipeline):
    condition_col = config["metadata"]["condition_column"]
    recovered = set(config["metadata"]["recovered_labels"])
    nonrecovered = set(config["metadata"]["nonrecovered_labels"])
    evaluate_conditions = set(config["perturbation"]["evaluate_conditions"])

    original_x = original_expr.loc[genes, meta.index].T
    perturbed_x = perturbed_expr.loc[genes, meta.index].T

    original_prob = pipeline.predict_proba(original_x)[:, 1]
    perturbed_prob = pipeline.predict_proba(perturbed_x)[:, 1]

    recovered_cells = meta.index[meta[condition_col].isin(recovered)]
    nonrecovered_cells = meta.index[meta[condition_col].isin(nonrecovered)]
    rec_centroid = original_x.loc[recovered_cells].mean(axis=0).values
    nonrec_centroid = original_x.loc[nonrecovered_cells].mean(axis=0).values

    orig_mat = original_x.values
    pert_mat = perturbed_x.values
    dist_orig_rec = distance_to_centroid(orig_mat, rec_centroid)
    dist_pert_rec = distance_to_centroid(pert_mat, rec_centroid)
    dist_orig_nonrec = distance_to_centroid(orig_mat, nonrec_centroid)
    dist_pert_nonrec = distance_to_centroid(pert_mat, nonrec_centroid)

    scores = pd.DataFrame({
        "model": model,
        "cell_id": meta.index,
        "condition": meta[condition_col].values,
        "original_recovery_probability": original_prob,
        "perturbed_recovery_probability": perturbed_prob,
        "delta_recovery_probability": perturbed_prob - original_prob,
        "dist_original_to_recovered_centroid": dist_orig_rec,
        "dist_perturbed_to_recovered_centroid": dist_pert_rec,
        "recovery_shift_score": dist_orig_rec - dist_pert_rec,
        "dist_original_to_nonrecovered_centroid": dist_orig_nonrec,
        "dist_perturbed_to_nonrecovered_centroid": dist_pert_nonrec,
        "nonrecovery_escape_score": dist_pert_nonrec - dist_orig_nonrec,
        "evaluated_post_cell": meta[condition_col].isin(evaluate_conditions).values,
    })

    delta = perturbed_expr.loc[genes, meta.index] - original_expr.loc[genes, meta.index]
    signature = pd.DataFrame({
        "model": model,
        "gene": genes,
        "mean_delta_expression": delta.mean(axis=1).values,
        "median_delta_expression": delta.median(axis=1).values,
        "sd_delta_expression": delta.std(axis=1).values,
    })
    return scores, signature


def summarize(scores):
    metrics = [
        "original_recovery_probability",
        "perturbed_recovery_probability",
        "delta_recovery_probability",
        "recovery_shift_score",
        "nonrecovery_escape_score",
    ]
    return (
        scores.groupby("model")[metrics].agg(["mean", "median", "std", "count"]),
        scores.groupby(["model", "condition"])[metrics].agg(["mean", "median", "std", "count"]),
        scores[scores["evaluated_post_cell"]].groupby("model")[metrics].agg(["mean", "median", "std", "count"]),
    )


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", default="configs/config.yaml")
    args = parser.parse_args()

    config = load_config(args.config)
    tables_dir = Path(config["paths"]["tables_dir"])
    tables_dir.mkdir(parents=True, exist_ok=True)
    prefix = config["model_results"].get("quantification_prefix", "model_quantification")

    manifest_path = Path(config["model_results"]["manifest_csv"])
    if not manifest_path.exists():
        raise SystemExit("Run scripts/01_check_and_transform_model_results.py first.")
    manifest = pd.read_csv(manifest_path)
    ready = manifest[manifest["status"] == "ready"]
    if ready.empty:
        print("No ready model results to quantify.")
        return

    expr = read_expression_matrix(config["paths"]["exported_expression_csv"])
    meta, expr = load_metadata(config, expr)
    genes = selected_genes(config, expr)
    pipeline = fit_classifier(config, meta, expr, genes)

    all_scores = []
    all_signatures = []
    for row in ready.to_dict("records"):
        output = read_expression_matrix(row["standard_path"])
        cells = meta.index.intersection(output.columns)
        if len(cells) == 0:
            print(f"Skipping {row['model']}: no overlapping cells")
            continue
        model_meta = meta.loc[cells].copy()
        perturbed = build_perturbed_matrix(output, expr, genes, cells, row["format"])
        scores, signature = score_model(
            config=config,
            model=row["model"],
            meta=model_meta,
            original_expr=expr,
            perturbed_expr=perturbed,
            genes=genes,
            pipeline=pipeline,
        )
        all_scores.append(scores)
        all_signatures.append(signature)

    if not all_scores:
        print("No model results could be quantified.")
        return

    scores = pd.concat(all_scores, ignore_index=True)
    signatures = pd.concat(all_signatures, ignore_index=True)
    scores.to_csv(tables_dir / f"{prefix}_cell_scores.csv", index=False)
    signatures.to_csv(tables_dir / f"{prefix}_gene_signatures.csv", index=False)

    overall, by_condition, post_only = summarize(scores)
    overall.to_csv(tables_dir / f"{prefix}_overall_summary.csv")
    by_condition.to_csv(tables_dir / f"{prefix}_by_condition_summary.csv")
    post_only.to_csv(tables_dir / f"{prefix}_post_only_summary.csv")

    print(f"Quantified {scores['model'].nunique()} model(s).")
    print(post_only.to_string())


if __name__ == "__main__":
    main()
