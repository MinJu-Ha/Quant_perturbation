import argparse
from pathlib import Path
import yaml
import numpy as np
import pandas as pd

from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import Pipeline
from sklearn.linear_model import LogisticRegression


def load_config(path):
    with open(path, "r") as f:
        return yaml.safe_load(f)


def read_expression_matrix(path):
    expr = pd.read_csv(path, index_col=0)
    expr.index = expr.index.astype(str)
    expr.columns = expr.columns.astype(str)
    return expr


def pick_target_gene(expr_index, candidates):
    for g in candidates:
        if g in expr_index:
            return g
    lower_map = {g.lower(): g for g in expr_index}
    for g in candidates:
        if g.lower() in lower_map:
            return lower_map[g.lower()]
    raise ValueError(f"Could not find target gene candidates: {candidates}")


def infer_target_signs(expr, target_gene, target_genes):
    # Infer whether each target is positively or negatively associated with BCL6
    # using cell-level Pearson correlation. This is a heuristic.
    bcl6 = expr.loc[target_gene].astype(float)
    signs = {}
    corrs = {}
    for g in target_genes:
        x = expr.loc[g].astype(float)
        corr = np.corrcoef(bcl6.values, x.values)[0, 1]
        if np.isnan(corr):
            corr = 0.0
        signs[g] = np.sign(corr)
        corrs[g] = corr
    return pd.DataFrame({
        "target_gene": list(signs.keys()),
        "correlation_with_BCL6": [corrs[g] for g in signs],
        "inferred_sign": [signs[g] for g in signs],
    }).set_index("target_gene")


def simulate_knockdown(expr, target_genes, signs_df, strength=0.5):
    # If a target is positively correlated with BCL6, knockdown decreases it.
    # If negatively correlated, knockdown increases it.
    perturbed = expr.copy()
    for g in target_genes:
        sign = signs_df.loc[g, "inferred_sign"]
        if sign == 0:
            continue
        sd = expr.loc[g].std()
        perturbed.loc[g] = expr.loc[g] - strength * sign * sd
    return perturbed


def distance_to_centroid(matrix_cells_by_genes, centroid):
    return np.linalg.norm(matrix_cells_by_genes - centroid[None, :], axis=1)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", default="configs/config.yaml")
    args = parser.parse_args()

    config = load_config(args.config)

    meta_path = Path(config["paths"]["exported_meta_csv"])
    expr_path = Path(config["paths"]["exported_expression_csv"])
    tables_dir = Path(config["paths"]["tables_dir"])
    tables_dir.mkdir(parents=True, exist_ok=True)

    meta = pd.read_csv(meta_path, index_col=0)
    expr = read_expression_matrix(expr_path)

    if "cell_id" in meta.columns:
        meta.index = meta["cell_id"].astype(str)

    common_cells = meta.index.astype(str).intersection(expr.columns.astype(str))
    meta = meta.loc[common_cells].copy()
    expr = expr.loc[:, common_cells].copy()

    condition_col = config["metadata"]["condition_column"]
    recovered = set(config["metadata"]["recovered_labels"])
    nonrecovered = set(config["metadata"]["nonrecovered_labels"])
    evaluate_conditions = set(config["perturbation"]["evaluate_conditions"])

    target_gene = pick_target_gene(expr.index, config["perturbation"]["target_gene_candidates"])

    target_file = tables_dir / "selected_bcl6_target_genes.csv"
    if not target_file.exists():
        raise FileNotFoundError("Run 02_bcl6_target_recovery_classifier.py first.")
    target_genes = pd.read_csv(target_file).iloc[:, 0].astype(str)
    target_genes = [g for g in target_genes if g in expr.index and g != target_gene]

    keep = meta[condition_col].isin(recovered.union(nonrecovered))
    meta_model = meta.loc[keep].copy()
    X = expr.loc[target_genes, meta_model.index].T
    y = meta_model[condition_col].isin(recovered).astype(int).values

    pipeline = Pipeline([
        ("scaler", StandardScaler(with_mean=True, with_std=True)),
        ("clf", LogisticRegression(max_iter=config["classifier"]["max_iter"], class_weight="balanced"))
    ])
    pipeline.fit(X, y)

    original_X_all = expr.loc[target_genes, meta.index].T
    original_prob = pipeline.predict_proba(original_X_all)[:, 1]

    signs_df = infer_target_signs(expr, target_gene, target_genes)
    signs_df.to_csv(tables_dir / "bcl6_target_inferred_signs.csv")

    perturbed_expr = simulate_knockdown(
        expr,
        target_genes=target_genes,
        signs_df=signs_df,
        strength=float(config["perturbation"]["knockdown_strength"]),
    )

    perturbed_X_all = perturbed_expr.loc[target_genes, meta.index].T
    perturbed_prob = pipeline.predict_proba(perturbed_X_all)[:, 1]

    recovered_cells = meta.index[meta[condition_col].isin(recovered)]
    nonrecovered_cells = meta.index[meta[condition_col].isin(nonrecovered)]

    rec_centroid = original_X_all.loc[recovered_cells].mean(axis=0).values
    nonrec_centroid = original_X_all.loc[nonrecovered_cells].mean(axis=0).values

    original_mat = original_X_all.values
    perturbed_mat = perturbed_X_all.values

    dist_orig_to_rec = distance_to_centroid(original_mat, rec_centroid)
    dist_pert_to_rec = distance_to_centroid(perturbed_mat, rec_centroid)

    dist_orig_to_nonrec = distance_to_centroid(original_mat, nonrec_centroid)
    dist_pert_to_nonrec = distance_to_centroid(perturbed_mat, nonrec_centroid)

    recovery_shift_score = dist_orig_to_rec - dist_pert_to_rec
    nonrecovery_escape_score = dist_pert_to_nonrec - dist_orig_to_nonrec

    out = pd.DataFrame({
        "cell_id": meta.index,
        "condition": meta[condition_col].values,
        "original_recovery_probability": original_prob,
        "perturbed_recovery_probability": perturbed_prob,
        "delta_recovery_probability": perturbed_prob - original_prob,
        "dist_original_to_recovered_centroid": dist_orig_to_rec,
        "dist_perturbed_to_recovered_centroid": dist_pert_to_rec,
        "recovery_shift_score": recovery_shift_score,
        "dist_original_to_nonrecovered_centroid": dist_orig_to_nonrec,
        "dist_perturbed_to_nonrecovered_centroid": dist_pert_to_nonrec,
        "nonrecovery_escape_score": nonrecovery_escape_score,
        "evaluated_post_cell": meta[condition_col].isin(evaluate_conditions).values,
    })

    out.to_csv(tables_dir / "bcl6_perturbation_cell_scores.csv", index=False)

    summary = out.groupby("condition")[[
        "original_recovery_probability",
        "perturbed_recovery_probability",
        "delta_recovery_probability",
        "recovery_shift_score",
        "nonrecovery_escape_score",
    ]].agg(["mean", "median", "std", "count"])

    summary.to_csv(tables_dir / "bcl6_perturbation_summary.csv")

    post = out[out["evaluated_post_cell"]]
    post_summary = post[[
        "delta_recovery_probability",
        "recovery_shift_score",
        "nonrecovery_escape_score"
    ]].agg(["mean", "median", "std", "count"])
    post_summary.to_csv(tables_dir / "bcl6_perturbation_post_only_summary.csv")

    print("Done.")
    print("Wrote results/tables/bcl6_perturbation_cell_scores.csv")


if __name__ == "__main__":
    main()
