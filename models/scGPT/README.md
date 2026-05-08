# scGPT

Notebook:

```text
run_scgpt_bcl6.ipynb
```

Checkpoint directory required by the notebook:

```text
SCGPT_CKPT_DIR/
  best_model.pt
  vocab.json
  args.json
```

Put one of these files in `results/`:

```text
scgpt_bcl6_knockdown_expression.csv
scgpt_bcl6_native.h5ad
scgpt_bcl6_native.npz
```

scGPT output is treated as `perturbed_expression`.

The `bcl6-scgpt` conda environment has already been created on this machine.
