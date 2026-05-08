# CellOracle

Notebook:

```text
run_celloracle_bcl6.ipynb
```

Put one of these files in `results/`:

```text
celloracle_bcl6_knockdown_delta.csv
celloracle_bcl6_native.h5ad
celloracle_bcl6_native.npz
```

CellOracle output is treated as `delta_expression`, so the evaluator uses:

```text
perturbed = original + delta
```

After the file exists, run:

```bash
python scripts/01_check_and_transform_model_results.py --config configs/config.yaml
```
