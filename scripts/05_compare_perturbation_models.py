import argparse
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd
import yaml
from sklearn.linear_model import LogisticRegression
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler


VALID_FORMATS = {"perturbed_expression", "delta_expression"}


def load_config(path):
    with open(path, "r") as f:
        return yaml.safe_load(f)


def read_expression_matrix(path):
    matrix = pd.read_csv(path, index_col=0)
    matrix.index = matrix.index.astype(str)
    matrix.columns = matrix.columns.astype(str)
    return matrix


def distance_to_centroid(matrix_cells_by_genes, centroid):
    return np.linalg.norm(matrix_cells_by_genes - centroid[None, :], axis=1)


def load_metadata(config, expr):
    meta = pd.read_csv(config["paths"]["exported_meta_csv"], index_col=0)
    if "cell_id" in meta.columns:
        meta.index = meta["cell_id"].astype(str)
    meta.index = meta.index.astype(str)

    common_cells = meta.index.intersection(expr.columns.astype(str))
    if len(common_cells) == 0:
        raise ValueError("No overlapping cell IDs between metadata and expression matrix.")
    return meta.loc[common_cells].copy(), expr.loc[:, common_cells].copy()


def load_selected_genes(config, expr):
    target_file = Path(config["paths"]["tables_dir"]) / "selected_bcl6_target_genes.csv"
    if not target_file.exists():
        raise FileNotFoundError("Run scripts/02_bcl6_target_recovery_classifier.py first.")

    genes = pd.read_csv(target_file).iloc[:, 0].dropna().astype(str)
    genes = [g for g in genes if g in expr.index]
    if len(genes) == 0:
        raise ValueError("No selected BCL6 target genes are present in the expression matrix.")
    return genes


def fit_recovery_classifier(config, meta, expr, genes):
    condition_col = config["metadata"]["condition_column"]
    recovered = set(config["metadata"]["recovered_labels"])
    nonrecovered = set(config["metadata"]["nonrecovered_labels"])

    keep = meta[condition_col].isin(recovered.union(nonrecovered))
    meta_model = meta.loc[keep].copy()
    X = expr.loc[genes, meta_model.index].T
    y = meta_model[condition_col].isin(recovered).astype(int).values

    if len(np.unique(y)) < 2:
        raise ValueError("Need both recovered and non-recovered labels for classifier.")

    pipeline = Pipeline([
        ("scaler", StandardScaler(with_mean=True, with_std=True)),
        ("clf", LogisticRegression(
            max_iter=config["classifier"]["max_iter"],
            class_weight="balanced",
        )),
    ])
    pipeline.fit(X, y)
    return pipeline


def prepare_perturbed_matrix(model_cfg, original_expr):
    path = Path(model_cfg["path"])
    if not path.exists():
        return None, f"missing file: {path}"

    output_format = model_cfg.get("format", "perturbed_expression")
    if output_format not in VALID_FORMATS:
        return None, f"invalid format {output_format!r}; expected one of {sorted(VALID_FORMATS)}"

    matrix = read_expression_matrix(path)
    common_genes = original_expr.index.intersection(matrix.index)
    common_cells = original_expr.columns.intersection(matrix.columns)
    if len(common_genes) == 0 or len(common_cells) == 0:
        return None, "no overlapping genes or cells with original expression matrix"

    aligned_original = original_expr.loc[common_genes, common_cells]
    aligned_output = matrix.loc[common_genes, common_cells]

    if output_format == "delta_expression":
        perturbed = aligned_original + aligned_output
    else:
        perturbed = aligned_output

    return perturbed, None


def score_model(config, model_name, meta, original_expr, perturbed_expr, genes, pipeline):
    condition_col = config["metadata"]["condition_column"]
    recovered = set(config["metadata"]["recovered_labels"])
    nonrecovered = set(config["metadata"]["nonrecovered_labels"])
    evaluate_conditions = set(config["perturbation"]["evaluate_conditions"])

    common_genes = pd.Index(genes).intersection(perturbed_expr.index)
    common_cells = meta.index.intersection(perturbed_expr.columns)
    if len(common_genes) == 0 or len(common_cells) == 0:
        raise ValueError(f"{model_name}: no overlapping selected genes or cells.")

    meta = meta.loc[common_cells].copy()
    original_X = original_expr.loc[common_genes, common_cells].T
    perturbed_X = perturbed_expr.loc[common_genes, common_cells].T

    original_prob = pipeline.predict_proba(original_X)[:, 1]
    perturbed_prob = pipeline.predict_proba(perturbed_X)[:, 1]

    recovered_cells = meta.index[meta[condition_col].isin(recovered)]
    nonrecovered_cells = meta.index[meta[condition_col].isin(nonrecovered)]
    rec_centroid = original_X.loc[recovered_cells].mean(axis=0).values
    nonrec_centroid = original_X.loc[nonrecovered_cells].mean(axis=0).values

    original_mat = original_X.values
    perturbed_mat = perturbed_X.values
    dist_orig_to_rec = distance_to_centroid(original_mat, rec_centroid)
    dist_pert_to_rec = distance_to_centroid(perturbed_mat, rec_centroid)
    dist_orig_to_nonrec = distance_to_centroid(original_mat, nonrec_centroid)
    dist_pert_to_nonrec = distance_to_centroid(perturbed_mat, nonrec_centroid)

    delta_expr = perturbed_expr.loc[common_genes, common_cells] - original_expr.loc[common_genes, common_cells]

    scores = pd.DataFrame({
        "model": model_name,
        "cell_id": meta.index,
        "condition": meta[condition_col].values,
        "original_recovery_probability": original_prob,
        "perturbed_recovery_probability": perturbed_prob,
        "delta_recovery_probability": perturbed_prob - original_prob,
        "dist_original_to_recovered_centroid": dist_orig_to_rec,
        "dist_perturbed_to_recovered_centroid": dist_pert_to_rec,
        "recovery_shift_score": dist_orig_to_rec - dist_pert_to_rec,
        "dist_original_to_nonrecovered_centroid": dist_orig_to_nonrec,
        "dist_perturbed_to_nonrecovered_centroid": dist_pert_to_nonrec,
        "nonrecovery_escape_score": dist_pert_to_nonrec - dist_orig_to_nonrec,
        "evaluated_post_cell": meta[condition_col].isin(evaluate_conditions).values,
    })

    signature = pd.DataFrame({
        "model": model_name,
        "gene": common_genes,
        "mean_delta_expression": delta_expr.mean(axis=1).values,
        "median_delta_expression": delta_expr.median(axis=1).values,
        "sd_delta_expression": delta_expr.std(axis=1).values,
    })
    return scores, signature


def summarize_scores(scores):
    metrics = [
        "original_recovery_probability",
        "perturbed_recovery_probability",
        "delta_recovery_probability",
        "recovery_shift_score",
        "nonrecovery_escape_score",
    ]
    by_condition = scores.groupby(["model", "condition"])[metrics].agg(["mean", "median", "std", "count"])
    post_only = scores[scores["evaluated_post_cell"]].groupby("model")[metrics].agg(["mean", "median", "std", "count"])
    overall = scores.groupby("model")[metrics].agg(["mean", "median", "std", "count"])
    return overall, by_condition, post_only


def pairwise_cell_concordance(scores):
    rows = []
    wide_delta = scores.pivot_table(
        index="cell_id",
        columns="model",
        values="delta_recovery_probability",
    )
    wide_shift = scores.pivot_table(
        index="cell_id",
        columns="model",
        values="recovery_shift_score",
    )

    for a, b in combinations(wide_delta.columns, 2):
        common_delta = wide_delta[[a, b]].dropna()
        common_shift = wide_shift[[a, b]].dropna()
        delta_corr = common_delta[a].corr(common_delta[b], method="spearman")
        shift_corr = common_shift[a].corr(common_shift[b], method="spearman")
        sign_agreement = (
            np.sign(common_delta[a].values) == np.sign(common_delta[b].values)
        ).mean()
        rows.append({
            "model_a": a,
            "model_b": b,
            "n_common_cells": len(common_delta),
            "spearman_delta_recovery_probability": delta_corr,
            "spearman_recovery_shift_score": shift_corr,
            "delta_sign_agreement": sign_agreement,
        })
    return pd.DataFrame(rows)


def pairwise_signature_concordance(signatures):
    rows = []
    wide = signatures.pivot_table(
        index="gene",
        columns="model",
        values="mean_delta_expression",
    )
    for a, b in combinations(wide.columns, 2):
        common = wide[[a, b]].dropna()
        rows.append({
            "model_a": a,
            "model_b": b,
            "n_common_genes": len(common),
            "pearson_mean_delta_expression": common[a].corr(common[b], method="pearson"),
            "spearman_mean_delta_expression": common[a].corr(common[b], method="spearman"),
            "delta_sign_agreement": (
                np.sign(common[a].values) == np.sign(common[b].values)
            ).mean(),
        })
    return pd.DataFrame(rows)


def write_input_template(config):
    out_dir = Path(config["paths"]["perturbation_outputs_dir"])
    out_dir.mkdir(parents=True, exist_ok=True)
    template = out_dir / "README_perturbation_output_format.md"
    template.write_text(
        "# Perturbation output format\n\n"
        "Provide one CSV per model. Rows are genes, columns are cells, and the first "
        "column contains gene names.\n\n"
        "Supported formats in `configs/config.yaml`:\n\n"
        "- `perturbed_expression`: values are the predicted post-perturbation expression.\n"
        "- `delta_expression`: values are predicted changes, so the evaluator uses "
        "`original_expression + delta_expression`.\n\n"
        "All values should be on the same scale as "
        "`data/processed/Cardiomyocyte_SCT_scaled.txt`.\n",
        encoding="utf-8",
    )
    return template


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", default="configs/config.yaml")
    parser.add_argument("--write-template", action="store_true")
    args = parser.parse_args()

    config = load_config(args.config)
    tables_dir = Path(config["paths"]["tables_dir"])
    tables_dir.mkdir(parents=True, exist_ok=True)

    if args.write_template:
        template = write_input_template(config)
        print(f"Wrote {template}")
        return

    expr = read_expression_matrix(config["paths"]["exported_expression_csv"])
    meta, expr = load_metadata(config, expr)
    genes = load_selected_genes(config, expr)
    pipeline = fit_recovery_classifier(config, meta, expr, genes)

    all_scores = []
    all_signatures = []
    skipped = []

    for model_cfg in config.get("model_comparison", {}).get("perturbation_outputs", []):
        model_name = model_cfg["name"]
        perturbed_expr, reason = prepare_perturbed_matrix(model_cfg, expr)
        if reason is not None:
            skipped.append({"model": model_name, "reason": reason})
            continue

        scores, signature = score_model(
            config=config,
            model_name=model_name,
            meta=meta,
            original_expr=expr,
            perturbed_expr=perturbed_expr,
            genes=genes,
            pipeline=pipeline,
        )
        all_scores.append(scores)
        all_signatures.append(signature)

    if skipped:
        pd.DataFrame(skipped).to_csv(tables_dir / "perturbation_model_skipped.csv", index=False)
        for item in skipped:
            print(f"Skipped {item['model']}: {item['reason']}")

    if not all_scores:
        print("No model perturbation outputs were available to compare.")
        print("Run with --write-template to create an input-format note.")
        return

    scores = pd.concat(all_scores, ignore_index=True)
    signatures = pd.concat(all_signatures, ignore_index=True)
    scores.to_csv(tables_dir / "perturbation_model_cell_scores.csv", index=False)
    signatures.to_csv(tables_dir / "perturbation_model_gene_signatures.csv", index=False)

    overall, by_condition, post_only = summarize_scores(scores)
    overall.to_csv(tables_dir / "perturbation_model_overall_summary.csv")
    by_condition.to_csv(tables_dir / "perturbation_model_by_condition_summary.csv")
    post_only.to_csv(tables_dir / "perturbation_model_post_only_summary.csv")

    pairwise_cell_concordance(scores).to_csv(
        tables_dir / "perturbation_model_pairwise_cell_concordance.csv",
        index=False,
    )
    pairwise_signature_concordance(signatures).to_csv(
        tables_dir / "perturbation_model_pairwise_gene_signature_concordance.csv",
        index=False,
    )

    print("Done.")
    print("Wrote cross-model perturbation comparison tables to results/tables.")


if __name__ == "__main__":
    main()
