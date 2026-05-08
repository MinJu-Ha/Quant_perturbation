# Geneformer

Notebook:

```text
run_geneformer_bcl6.ipynb
```

Prepared inputs:

```text
data/processed/symbol_to_ensembl.csv
data/processed/perturbations/cardiomyocyte.loom
```

Geneformer package/model still needs to be available locally. The model folder
should be passed via `GENEFORMER_DIR`.

Put one of these files in `results/`:

```text
geneformer_bcl6_knockdown_expression.csv
geneformer_bcl6_native.npz
```

`geneformer_bcl6_native.npz` must contain:

```text
cell_ids
original_embeddings
perturbed_embeddings
```

Geneformer output is treated as `perturbed_expression`.
