import argparse
import subprocess
from pathlib import Path

import yaml


def load_config(path):
    with open(path, "r") as f:
        return yaml.safe_load(f)


def run_step(name, command, required=True):
    print(f"\n===== {name} =====")
    print(" ".join(command))
    completed = subprocess.run(command, check=False)
    if completed.returncode != 0:
        message = f"{name} failed with exit code {completed.returncode}."
        if required:
            raise RuntimeError(message)
        print(message)
    return completed.returncode


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", default="configs/config.yaml")
    parser.add_argument("--skip-models", action="store_true")
    parser.add_argument("--skip-inspection", action="store_true")
    args = parser.parse_args()

    config = load_config(args.config)
    e2e = config.get("end_to_end", {})
    config_path = args.config

    Path(config["paths"]["processed_dir"]).mkdir(parents=True, exist_ok=True)
    Path(config["paths"]["tables_dir"]).mkdir(parents=True, exist_ok=True)
    Path(config["paths"]["figures_dir"]).mkdir(parents=True, exist_ok=True)
    Path(config["paths"]["perturbation_outputs_dir"]).mkdir(parents=True, exist_ok=True)

    if e2e.get("run_inspection", True) and not args.skip_inspection:
        run_step("Inspect Seurat object", ["Rscript", "scripts/00_inspect_post_cardiomyocyte.R"])

    if e2e.get("run_export", True):
        run_step("Export Seurat data", ["Rscript", "scripts/01_export_from_seurat.R"])

    if e2e.get("run_classifier", True):
        run_step(
            "Train recovery classifier",
            ["python", "scripts/02_bcl6_target_recovery_classifier.py", "--config", config_path],
        )

    if e2e.get("run_model_perturbations", True) and not args.skip_models:
        run_step(
            "Run model perturbations",
            ["python", "scripts/06_run_model_perturbations.py", "--config", config_path],
            required=e2e.get("stop_on_model_failure", False),
        )

    if e2e.get("run_comparison", True):
        run_step(
            "Compare perturbation outputs",
            ["python", "scripts/05_compare_perturbation_models.py", "--config", config_path],
            required=False,
        )

    print("\nPipeline finished.")


if __name__ == "__main__":
    main()
