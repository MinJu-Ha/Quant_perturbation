import argparse
from pathlib import Path
import yaml
import pandas as pd
import matplotlib.pyplot as plt


def load_config(path):
    with open(path, "r") as f:
        return yaml.safe_load(f)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", default="configs/config.yaml")
    args = parser.parse_args()

    config = load_config(args.config)

    tables_dir = Path(config["paths"]["tables_dir"])
    figures_dir = Path(config["paths"]["figures_dir"])
    figures_dir.mkdir(parents=True, exist_ok=True)

    scores_path = tables_dir / "bcl6_perturbation_cell_scores.csv"
    if not scores_path.exists():
        raise FileNotFoundError("Run 03_simulate_bcl6_perturbation.py first.")

    df = pd.read_csv(scores_path)

    plt.figure(figsize=(7, 5))
    plt.hist(df["delta_recovery_probability"].dropna(), bins=60)
    plt.axvline(0, linestyle="--")
    plt.xlabel("Delta recovery probability")
    plt.ylabel("Number of cells")
    plt.title("BCL6 perturbation: recovery probability shift")
    plt.tight_layout()
    plt.savefig(figures_dir / "recovery_probability_shift.png", dpi=300)
    plt.close()

    conditions = list(df["condition"].dropna().unique())
    data = [df.loc[df["condition"] == c, "delta_recovery_probability"].dropna() for c in conditions]

    plt.figure(figsize=(8, 5))
    plt.boxplot(data, labels=conditions, showfliers=False)
    plt.axhline(0, linestyle="--")
    plt.ylabel("Delta recovery probability")
    plt.title("BCL6 perturbation shift by condition")
    plt.xticks(rotation=30, ha="right")
    plt.tight_layout()
    plt.savefig(figures_dir / "recovery_shift_by_condition.png", dpi=300)
    plt.close()

    data2 = [df.loc[df["condition"] == c, "recovery_shift_score"].dropna() for c in conditions]

    plt.figure(figsize=(8, 5))
    plt.boxplot(data2, labels=conditions, showfliers=False)
    plt.axhline(0, linestyle="--")
    plt.ylabel("Recovery centroid shift score")
    plt.title("BCL6 perturbation: distance shift toward recovered centroid")
    plt.xticks(rotation=30, ha="right")
    plt.tight_layout()
    plt.savefig(figures_dir / "recovery_centroid_shift_by_condition.png", dpi=300)
    plt.close()

    print("Done. Check results/figures.")


if __name__ == "__main__":
    main()
