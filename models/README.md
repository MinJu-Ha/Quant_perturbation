# Model Perturbation Runs

Each model owns its own folder:

- `CellOracle/`
- `scGPT/`
- `Geneformer/`

Run or edit the notebook inside each folder, then put that model's output in
its local `results/` directory. The project-level scripts only read from these
folders and skip any model whose result is missing.

The KO gene is configured in `configs/config.yaml` under `perturbation.gene`.

Expected standard output:

```text
models/<Model>/results/standard_perturbation.csv
```

Rows are genes, columns are cells. Values must be on the same scale as
`paths.exported_expression_csv`.

Pipeline:

```bash
python scripts/01_check_and_transform_model_results.py --config configs/config.yaml
python scripts/02_quantify_model_results.py --config configs/config.yaml
python scripts/03_compare_and_plot_model_results.py --config configs/config.yaml
```
