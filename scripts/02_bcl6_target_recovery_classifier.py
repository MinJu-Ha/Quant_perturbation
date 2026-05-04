import argparse
from pathlib import Path
import yaml
import numpy as np
import pandas as pd

from sklearn.model_selection import train_test_split, StratifiedKFold, cross_val_score
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import Pipeline
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score, average_precision_score, accuracy_score


def load_config(path):
    with open(path, "r") as f:
        return yaml.safe_load(f)


def read_expression_matrix(path):
    expr = pd.read_csv(path, index_col=0)
    expr.index = expr.index.astype(str)
    expr.columns = expr.columns.astype(str)
    return expr


def read_chip_targets(path):
    chip = pd.read_csv(path)
    first_col = chip.columns[0]
    genes = chip[first_col].dropna().astype(str).str.replace("\ufeff", "", regex=False)
    return pd.Index(genes.unique())


def pick_target_gene(expr_index, candidates):
    for g in candidates:
        if g in expr_index:
            return g
    lower_map = {g.lower(): g for g in expr_index}
    for g in candidates:
        if g.lower() in lower_map:
            return lower_map[g.lower()]
    raise ValueError(f"Could not find target gene candidates in expression matrix: {candidates}")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", default="configs/config.yaml")
    args = parser.parse_args()

    config = load_config(args.config)

    meta_path = Path(config["paths"]["exported_meta_csv"])
    expr_path = Path(config["paths"]["exported_expression_csv"])
    chip_path = Path(config["paths"]["chip_targets_csv"])
    tables_dir = Path(config["paths"]["tables_dir"])
    tables_dir.mkdir(parents=True, exist_ok=True)

    meta = pd.read_csv(meta_path, index_col=0)
    expr = read_expression_matrix(expr_path)
    chip_targets = read_chip_targets(chip_path)

    condition_col = config["metadata"]["condition_column"]
    recovered = set(config["metadata"]["recovered_labels"])
    nonrecovered = set(config["metadata"]["nonrecovered_labels"])

    common_cells = meta.index.astype(str).intersection(expr.columns.astype(str))
    if len(common_cells) == 0 and "cell_id" in meta.columns:
        meta = meta.set_index(meta["cell_id"].astype(str), drop=False)
        common_cells = meta.index.astype(str).intersection(expr.columns.astype(str))
    if len(common_cells) == 0:
        raise ValueError("No overlapping cell IDs between metadata and expression matrix.")

    meta = meta.loc[common_cells].copy()
    expr = expr.loc[:, common_cells].copy()

    target_gene = pick_target_gene(expr.index, config["perturbation"]["target_gene_candidates"])

    exact_targets = pd.Index([g for g in chip_targets if g in expr.index])
    expr_upper_to_gene = {g.upper(): g for g in expr.index}
    ci_targets = []
    for g in chip_targets:
        if g.upper() in expr_upper_to_gene:
            ci_targets.append(expr_upper_to_gene[g.upper()])
    ci_targets = pd.Index(pd.unique(ci_targets))

    targets = exact_targets if len(exact_targets) >= config["features"]["min_target_genes"] else ci_targets
    targets = pd.Index([g for g in targets if g != target_gene])

    if len(targets) < config["features"]["min_target_genes"]:
        raise ValueError(f"Too few BCL6 targets found in expression matrix: {len(targets)}")

    max_targets = int(config["features"]["max_target_genes"])
    if len(targets) > max_targets:
        variances = expr.loc[targets].var(axis=1).sort_values(ascending=False)
        targets = pd.Index(variances.head(max_targets).index)

    pd.DataFrame({"target_gene": targets}).to_csv(tables_dir / "bcl6_target_overlap.csv", index=False)

    keep = meta[condition_col].isin(recovered.union(nonrecovered))
    meta_model = meta.loc[keep].copy()
    X = expr.loc[targets, meta_model.index].T
    y = meta_model[condition_col].isin(recovered).astype(int).values

    if len(np.unique(y)) < 2:
        raise ValueError("Need both recovered and non-recovered labels for classifier.")

    pipeline = Pipeline([
        ("scaler", StandardScaler(with_mean=True, with_std=True)),
        ("clf", LogisticRegression(max_iter=config["classifier"]["max_iter"], class_weight="balanced"))
    ])

    X_train, X_test, y_train, y_test = train_test_split(
        X, y,
        test_size=config["classifier"]["test_size"],
        random_state=config["classifier"]["random_state"],
        stratify=y
    )

    pipeline.fit(X_train, y_train)
    prob_test = pipeline.predict_proba(X_test)[:, 1]
    pred_test = (prob_test >= 0.5).astype(int)

    metrics = {
        "n_cells_total": len(meta_model),
        "n_features_bcl6_targets": len(targets),
        "target_gene_used": target_gene,
        "test_roc_auc": roc_auc_score(y_test, prob_test),
        "test_average_precision": average_precision_score(y_test, prob_test),
        "test_accuracy": accuracy_score(y_test, pred_test),
    }

    cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=config["classifier"]["random_state"])
    cv_auc = cross_val_score(pipeline, X, y, cv=cv, scoring="roc_auc")
    metrics["cv_auc_mean"] = float(cv_auc.mean())
    metrics["cv_auc_std"] = float(cv_auc.std())

    pd.DataFrame([metrics]).to_csv(tables_dir / "recovery_classifier_metrics.csv", index=False)

    pipeline.fit(X, y)
    X_all = expr.loc[targets, meta.index].T
    recovery_prob = pipeline.predict_proba(X_all)[:, 1]

    score_df = pd.DataFrame({
        "cell_id": meta.index,
        "condition": meta[condition_col].values,
        "original_recovery_probability": recovery_prob,
        "is_labeled_recovered": meta[condition_col].isin(recovered).values,
        "is_labeled_nonrecovered": meta[condition_col].isin(nonrecovered).values,
    })
    score_df.to_csv(tables_dir / "original_recovery_probabilities.csv", index=False)

    pd.Series(targets, name="selected_bcl6_target_gene").to_csv(
        tables_dir / "selected_bcl6_target_genes.csv",
        index=False
    )

    print("Done.")
    print(metrics)


if __name__ == "__main__":
    main()
