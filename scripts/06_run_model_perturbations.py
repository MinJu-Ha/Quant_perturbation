import argparse
from pathlib import Path

import pandas as pd
import yaml

from perturbation_adapters.celloracle_adapter import CellOracleAdapter
from perturbation_adapters.external import ExternalCommandAdapter
from perturbation_adapters.geneformer_adapter import GeneformerAdapter
from perturbation_adapters.scgpt_adapter import ScGPTAdapter


PYTHON_API_ADAPTERS = {
    "celloracle": CellOracleAdapter,
    "scgpt": ScGPTAdapter,
    "geneformer": GeneformerAdapter,
}


def load_config(path):
    with open(path, "r") as f:
        return yaml.safe_load(f)


def run_model(config, model_name, runner_cfg):
    mode = runner_cfg.get("mode", "external_command")
    if mode == "external_command":
        adapter = ExternalCommandAdapter(config, runner_cfg)
        adapter.model_name = model_name
    elif mode == "python_api":
        adapter_cls = PYTHON_API_ADAPTERS.get(model_name)
        if adapter_cls is None:
            raise ValueError(f"No python_api adapter registered for {model_name}.")
        adapter = adapter_cls(config, runner_cfg)
    else:
        raise ValueError(f"Unknown mode for {model_name}: {mode}")

    return adapter.run()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", default="configs/config.yaml")
    parser.add_argument("--models", nargs="*", default=None, help="Optional subset: celloracle scgpt geneformer")
    args = parser.parse_args()

    config = load_config(args.config)
    runners = config.get("model_runners", {})
    selected = set(args.models) if args.models else set(runners)
    tables_dir = Path(config["paths"]["tables_dir"])
    tables_dir.mkdir(parents=True, exist_ok=True)

    results = []
    for model_name, runner_cfg in runners.items():
        if model_name not in selected:
            continue
        if not runner_cfg.get("enabled", False):
            results.append({
                "model": model_name,
                "status": "skipped",
                "message": "disabled in config",
                "output_path": runner_cfg.get("output_path"),
                "output_format": runner_cfg.get("output_format"),
            })
            continue

        result = run_model(config, model_name, runner_cfg)
        results.append({
            "model": result.model,
            "status": result.status,
            "message": result.message,
            "output_path": str(result.output_path),
            "output_format": result.output_format,
        })

        if result.status != "completed" and config.get("end_to_end", {}).get("stop_on_model_failure", False):
            break

    out = pd.DataFrame(results)
    out.to_csv(tables_dir / "perturbation_model_run_status.csv", index=False)
    print(out.to_string(index=False))


if __name__ == "__main__":
    main()
