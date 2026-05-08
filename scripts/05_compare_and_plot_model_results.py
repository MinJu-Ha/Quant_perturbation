import argparse
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd
import yaml


def load_config(path):
    with open(path, "r") as f:
        return yaml.safe_load(f)


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
        rows.append({
            "model_a": a,
            "model_b": b,
            "n_common_cells": len(common_delta),
            "spearman_delta_recovery_probability": common_delta[a].corr(common_delta[b], method="spearman"),
            "spearman_recovery_shift_score": common_shift[a].corr(common_shift[b], method="spearman"),
            "delta_sign_agreement": (
                np.sign(common_delta[a].values) == np.sign(common_delta[b].values)
            ).mean() if len(common_delta) else np.nan,
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
            ).mean() if len(common) else np.nan,
        })
    return pd.DataFrame(rows)


def plot_box(df, value, by, output, title, ylabel):
    import matplotlib.pyplot as plt

    groups = list(df[by].dropna().unique())
    data = [df.loc[df[by] == group, value].dropna() for group in groups]
    plt.figure(figsize=(max(6, len(groups) * 1.4), 5))
    plt.boxplot(data, tick_labels=groups, showfliers=False)
    plt.axhline(0, linestyle="--", color="black", linewidth=1)
    plt.title(title)
    plt.ylabel(ylabel)
    plt.xticks(rotation=30, ha="right")
    plt.tight_layout()
    plt.savefig(output, dpi=300)
    plt.close()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", default="configs/config.yaml")
    args = parser.parse_args()

    config = load_config(args.config)
    tables_dir = Path(config["paths"]["tables_dir"])
    figures_dir = Path(config["paths"]["figures_dir"])
    figures_dir.mkdir(parents=True, exist_ok=True)
    prefix = config["model_results"].get("quantification_prefix", "model_quantification")

    scores_path = tables_dir / f"{prefix}_cell_scores.csv"
    signatures_path = tables_dir / f"{prefix}_gene_signatures.csv"
    if not scores_path.exists():
        print("No quantified model results found. Run scripts/02_quantify_model_results.py after at least one model result is ready.")
        return

    scores = pd.read_csv(scores_path)
    signatures = pd.read_csv(signatures_path) if signatures_path.exists() else pd.DataFrame()

    pairwise_cells = pairwise_cell_concordance(scores)
    pairwise_cells.to_csv(tables_dir / f"{prefix}_pairwise_cell_concordance.csv", index=False)

    if not signatures.empty:
        pairwise_signatures = pairwise_signature_concordance(signatures)
    else:
        pairwise_signatures = pd.DataFrame()
    pairwise_signatures.to_csv(tables_dir / f"{prefix}_pairwise_gene_signature_concordance.csv", index=False)

    plot_box(
        scores,
        value="delta_recovery_probability",
        by="model",
        output=figures_dir / f"{prefix}_delta_recovery_by_model.png",
        title="KO perturbation recovery probability shift by model",
        ylabel="Delta recovery probability",
    )

    post = scores[scores["evaluated_post_cell"]].copy()
    if not post.empty:
        plot_box(
            post,
            value="delta_recovery_probability",
            by="model",
            output=figures_dir / f"{prefix}_post_delta_recovery_by_model.png",
            title="Post-cell KO perturbation shift by model",
            ylabel="Delta recovery probability",
        )

    plot_box(
        scores,
        value="recovery_shift_score",
        by="model",
        output=figures_dir / f"{prefix}_recovery_centroid_shift_by_model.png",
        title="KO perturbation centroid shift by model",
        ylabel="Recovery centroid shift score",
    )

    print(f"Wrote comparison tables and plots for {scores['model'].nunique()} model(s).")


if __name__ == "__main__":
    main()
