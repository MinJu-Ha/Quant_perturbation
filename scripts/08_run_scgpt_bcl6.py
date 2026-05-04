"""Convert an scGPT native perturbation result into the standard
genes x cells perturbed-expression CSV consumed by
scripts/05_compare_perturbation_models.py.

Two ways to use this script:

1. End-to-end: run notebooks/20_scgpt_bcl6.ipynb to produce the native
   artefact, then run this script.

2. Bring-your-own-result: if you already have an AnnData (.h5ad) or .npz with
   the perturbed expression matrix from your own scGPT pipeline, point
   `model_runners.scgpt.native_output_path` at it (or pass --native-output)
   and run this script directly. No notebook needed.
"""
import argparse
import sys
from pathlib import Path

import yaml

sys.path.insert(0, str(Path(__file__).resolve().parent))
from converters import scgpt_native  # noqa: E402


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
            "Path to an scGPT native artefact (.h5ad/.npz/.csv). "
            "Falls back to model_runners.scgpt.native_output_path."
        ),
    )
    args = parser.parse_args()

    config = load_config(args.config)
    runner = config["model_runners"]["scgpt"]
    output_path = Path(runner["output_path"])
    output_path.parent.mkdir(parents=True, exist_ok=True)

    native_path = args.native_output or runner.get("native_output_path")
    if not native_path:
        raise SystemExit(
            "No scGPT native output configured. Either pass --native-output PATH "
            "or set model_runners.scgpt.native_output_path in configs/config.yaml. "
            "Run notebooks/20_scgpt_bcl6.ipynb to produce one."
        )
    native_path = Path(native_path)
    if not native_path.exists():
        raise SystemExit(
            f"scGPT native output not found at {native_path}. "
            "Run notebooks/20_scgpt_bcl6.ipynb first, or fix the configured path."
        )

    n_genes, n_cells = scgpt_native.convert(native_path, output_path)
    print(
        f"[scgpt] converted {native_path}\n"
        f"        -> {output_path} ({n_genes} genes x {n_cells} cells, perturbed_expression)"
    )


if __name__ == "__main__":
    main()
