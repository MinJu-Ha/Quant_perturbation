import argparse
import shutil
import sys
from pathlib import Path

import pandas as pd
import yaml


SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_DIR))

from converters import celloracle_native, geneformer_native, scgpt_native  # noqa: E402


CONVERTERS = {
    "CellOracle": celloracle_native.convert,
    "scGPT": scgpt_native.convert,
}


def load_config(path):
    with open(path, "r") as f:
        return yaml.safe_load(f)


def ko_gene(config):
    gene = config.get("perturbation", {}).get("gene")
    if gene:
        return str(gene)
    return str(config["perturbation"]["target_gene_candidates"][0])


def read_csv_axes(path):
    header = pd.read_csv(path, nrows=0)
    cells = pd.Index(header.columns[1:].astype(str))
    genes = pd.read_csv(path, usecols=[0]).iloc[:, 0].dropna().astype(str)
    return pd.Index(genes), cells


def pick_existing(folder, candidates):
    for rel in candidates:
        path = Path(folder) / rel
        if path.exists():
            return path
    return None


def convert_native(model_name, native_path, standard_path, config):
    suffix = native_path.suffix.lower()
    if suffix in {".csv", ".tsv", ".txt"}:
        standard_path.parent.mkdir(parents=True, exist_ok=True)
        if native_path.resolve() != standard_path.resolve():
            shutil.copyfile(native_path, standard_path)
        return

    if model_name == "Geneformer":
        ridge_alpha = 1.0
        for entry in config.get("model_results", {}).get("models", []):
            if entry.get("name") == "Geneformer":
                ridge_alpha = float(entry.get("ridge_alpha", 1.0))
                break
        geneformer_native.convert(
            native_path=native_path,
            output_path=standard_path,
            original_expression_csv=config["paths"]["exported_expression_csv"],
            ridge_alpha=ridge_alpha,
        )
        return

    converter = CONVERTERS.get(model_name)
    if converter is None:
        raise ValueError(f"No converter registered for {model_name}.")
    converter(native_path, standard_path)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", default="configs/config.yaml")
    args = parser.parse_args()

    config = load_config(args.config)
    model_cfg = config["model_results"]
    tables_dir = Path(config["paths"]["tables_dir"])
    tables_dir.mkdir(parents=True, exist_ok=True)

    gene = ko_gene(config)
    original_genes, original_cells = read_csv_axes(config["paths"]["exported_expression_csv"])
    gene_upper = gene.upper()
    original_gene_upper = {g.upper() for g in original_genes}
    ko_in_original = gene_upper in original_gene_upper

    rows = []
    for item in model_cfg["models"]:
        name = item["name"]
        folder = Path(item["folder"])
        result_dir = folder / "results"
        standard_path = result_dir / model_cfg.get("standard_filename", "standard_perturbation.csv")
        native_path = pick_existing(folder, item.get("native_candidates", []))

        row = {
            "model": name,
            "status": "missing",
            "reason": "",
            "native_path": str(native_path) if native_path else "",
            "standard_path": str(standard_path),
            "format": item["format"],
            "ko_gene": gene,
            "ko_gene_in_original_expression": ko_in_original,
            "ko_gene_in_model_output": False,
            "n_genes": 0,
            "n_cells": 0,
            "n_gene_overlap_with_original": 0,
            "n_cell_overlap_with_original": 0,
        }

        if native_path is None and not standard_path.exists():
            row["reason"] = f"no result file found in {result_dir}"
            rows.append(row)
            continue

        try:
            if native_path is not None:
                convert_native(name, native_path, standard_path, config)
            genes, cells = read_csv_axes(standard_path)
            gene_upper_set = {g.upper() for g in genes}
            row.update({
                "status": "ready",
                "reason": "",
                "ko_gene_in_model_output": gene_upper in gene_upper_set,
                "n_genes": len(genes),
                "n_cells": len(cells),
                "n_gene_overlap_with_original": len(original_genes.intersection(genes)),
                "n_cell_overlap_with_original": len(original_cells.intersection(cells)),
            })
            if row["n_gene_overlap_with_original"] == 0 or row["n_cell_overlap_with_original"] == 0:
                row["status"] = "invalid"
                row["reason"] = "no overlapping genes or cells with exported expression"
        except Exception as exc:
            row["status"] = "error"
            row["reason"] = str(exc)
        rows.append(row)

    manifest = pd.DataFrame(rows)
    manifest_path = Path(model_cfg["manifest_csv"])
    manifest_path.parent.mkdir(parents=True, exist_ok=True)
    manifest.to_csv(manifest_path, index=False)
    print(manifest.to_string(index=False))
    print(f"\nWrote {manifest_path}")


if __name__ == "__main__":
    main()
