"""Convert a Geneformer native perturbation result into the standard
genes x cells perturbed-expression CSV consumed by
scripts/05_compare_perturbation_models.py.

Geneformer natively returns embedding shifts. This script projects them back
into expression space using a ridge-regression linear map fitted on the
exported original expression matrix.

Two ways to use this script:

1. End-to-end: run notebooks/30_geneformer_bcl6.ipynb to produce the native
   .npz of original/perturbed embeddings, then run this script.

2. Bring-your-own-result: if you already have Geneformer original and perturbed
   embeddings from your own pipeline, save them as a .npz with three arrays
   ("cell_ids", "original_embeddings", "perturbed_embeddings"), point
   `model_runners.geneformer.native_output_path` at it (or pass
   --native-output) and run this script directly.
"""
import argparse
import sys
from pathlib import Path

import yaml

sys.path.insert(0, str(Path(__file__).resolve().parent))
from converters import geneformer_native  # noqa: E402


def load_config(path):
    with open(path, "r") as f:
        return yaml.safe_load(f)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", default="configs/config.yaml")
    parser.add_argument(
        "--native-output",
        default=None,
        help=(
            "Path to a Geneformer native embedding archive (.npz). "
            "Falls back to model_runners.geneformer.native_output_path."
        ),
    )
    parser.add_argument(
        "--ridge-alpha",
        type=float,
        default=None,
        help=(
            "Ridge regularisation strength for the embedding -> expression "
            "projection. Default 1.0, can also be set via "
            "model_runners.geneformer.ridge_alpha."
        ),
    )
    args = parser.parse_args()

    config = load_config(args.config)
    runner = config["model_runners"]["geneformer"]
    output_path = Path(runner["output_path"])
    output_path.parent.mkdir(parents=True, exist_ok=True)

    native_path = args.native_output or runner.get("native_output_path")
    if not native_path:
        raise SystemExit(
            "No Geneformer native output configured. Either pass --native-output PATH "
            "or set model_runners.geneformer.native_output_path in configs/config.yaml."
        )
    native_path = Path(native_path)
    if not native_path.exists():
        raise SystemExit(
            f"Geneformer native output not found at {native_path}. "
            "Run notebooks/30_geneformer_bcl6.ipynb first, or fix the configured path."
        )

    original_expression_csv = config["paths"]["exported_expression_csv"]
    if not Path(original_expression_csv).exists():
        raise SystemExit(
            f"Exported expression CSV not found at {original_expression_csv}. "
            "Run scripts/01_export_from_seurat.R first."
        )

    ridge_alpha = (
        args.ridge_alpha
        if args.ridge_alpha is not None
        else float(runner.get("ridge_alpha", 1.0))
    )

    n_genes, n_cells = geneformer_native.convert(
        native_path=native_path,
        output_path=output_path,
        original_expression_csv=original_expression_csv,
        ridge_alpha=ridge_alpha,
    )
    print(
        f"[geneformer] converted {native_path}\n"
        f"             -> {output_path} ({n_genes} genes x {n_cells} cells, "
        f"perturbed_expression, ridge_alpha={ridge_alpha})"
    )


if __name__ == "__main__":
    main()
