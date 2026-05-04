# BCL6 Recovery-State Perturbation Project

## Goal

Quantitatively test whether **BCL6 in silico perturbation** shifts post-LVAD cardiomyocytes toward a **recovery-like transcriptional state**.

The project is also set up to compare perturbation outputs from multiple models,
especially **CellOracle**, **scGPT**, and **Geneformer**, using the same recovery
classifier and centroid-distance metrics.

## Main question

For `post_cardiomyocyte.rds`:

> If we perturb BCL6 in post-LVAD cardiomyocytes, do cells become more similar to recovered cardiomyocytes than non-recovered cardiomyocytes?

## Data expected

Place these files in `data/raw/`:

```text
post_cardiomyocyte.rds
Bcl6_ChIP_target_P0.05.csv
```

Optional, if already exported:

```text
Cardiomyocyte_meta.csv
Cardiomyocyte_SCT_scaled.txt
```

## Labels

```text
Recovered-like: Rpre, Rpost
Non-recovered-like: NRpre, NRpost
Reference/control: Donor
```

## Quantitative outputs

1. Recovery classifier probability
2. Recovery signature score / centroid distance
3. BCL6 perturbation-induced shift
4. Cross-model concordance across CellOracle, scGPT, and Geneformer

Main score:

```text
delta_recovery_probability =
perturbed_recovery_probability - original_recovery_probability
```

Interpretation:

```text
delta > 0: BCL6 perturbation makes cells more recovery-like
delta < 0: BCL6 perturbation makes cells less recovery-like
```

Centroid score:

```text
recovery_shift_score =
distance(original cell, recovered centroid)
-
distance(perturbed cell, recovered centroid)
```

Interpretation:

```text
score > 0: perturbation moves the cell closer to recovered state
score < 0: perturbation moves the cell farther from recovered state
```

## Quick start

```bash
Rscript scripts/00_inspect_post_cardiomyocyte.R
Rscript scripts/01_export_from_seurat.R
python scripts/02_bcl6_target_recovery_classifier.py --config configs/config.yaml
python scripts/03_simulate_bcl6_perturbation.py --config configs/config.yaml
python scripts/04_plot_bcl6_perturbation_results.py --config configs/config.yaml
```

## Comparing CellOracle, scGPT, and Geneformer

This repo treats the three model runs as upstream perturbation engines. After
each model predicts a BCL6 perturbation result, export one matrix per model into
`data/processed/perturbations/`.

Expected format:

```text
gene,cell_1,cell_2,cell_3
GeneA,0.12,0.08,-0.03
GeneB,-0.31,-0.10,0.22
```

Rows are genes, columns are cells, and values must use the same expression scale
as `data/processed/Cardiomyocyte_SCT_scaled.txt`.

Supported model-output types in `configs/config.yaml`:

```text
perturbed_expression: predicted expression after BCL6 perturbation
delta_expression: predicted expression change after BCL6 perturbation
```

Default expected files:

```text
data/processed/perturbations/celloracle_bcl6_knockdown_delta.csv
data/processed/perturbations/scgpt_bcl6_knockdown_expression.csv
data/processed/perturbations/geneformer_bcl6_knockdown_expression.csv
```

Run:

```bash
python scripts/05_compare_perturbation_models.py --config configs/config.yaml --write-template
python scripts/05_compare_perturbation_models.py --config configs/config.yaml
```

Main comparison outputs:

```text
results/tables/perturbation_model_cell_scores.csv
results/tables/perturbation_model_overall_summary.csv
results/tables/perturbation_model_by_condition_summary.csv
results/tables/perturbation_model_post_only_summary.csv
results/tables/perturbation_model_pairwise_cell_concordance.csv
results/tables/perturbation_model_pairwise_gene_signature_concordance.csv
```

The key numerical comparison is still `delta_recovery_probability`, but now it is
computed separately for each model. The pairwise concordance tables quantify
whether CellOracle, scGPT, and Geneformer agree at the cell level and gene
perturbation-signature level.

## End-to-end workflow

The full intended workflow is:

```text
RDS + BCL6 target list
  -> inspect and export Seurat data
  -> train a recovered vs non-recovered classifier
  -> run the same BCL6 deletion/knockdown with CellOracle, scGPT, Geneformer
  -> collect each model output
  -> score each model with the same recovery metrics
  -> compare model-level and pairwise concordance results
```

Run the orchestrator:

```bash
python scripts/10_run_end_to_end_pipeline.py --config configs/config.yaml
```

For a dry run of the shared preprocessing/scoring path without running model
engines:

```bash
python scripts/10_run_end_to_end_pipeline.py --config configs/config.yaml --skip-models
```

Model execution is controlled by `model_runners` in `configs/config.yaml`.
Each runner has a **native artefact path** (the model's own raw output) and a
**standard output path** (the project's genes x cells CSV). The runner scripts
do nothing except convert from one to the other:

```text
notebook (or your own pipeline)
  -> writes <native_output_path>
  -> scripts/07-09 read it and write <output_path>
  -> scripts/05_compare_perturbation_models.py reads <output_path>
```

This split means a user who already has a CellOracle Oracle, an scGPT
perturbation `.h5ad`, or Geneformer embedding pickles only has to drop the
file at `native_output_path` and run scripts/07-09. They do not need to use
the notebooks.

### Per-model setup

Each model has its own conda environment so dependencies do not collide:

```bash
conda env create -f envs/celloracle.yml
conda env create -f envs/scgpt.yml
conda env create -f envs/geneformer.yml
```

Notebooks (one per model) reproduce the official tutorials, applied to this
project's exported data:

```text
notebooks/10_celloracle_bcl6.ipynb   -> data/processed/perturbations/celloracle_bcl6_native.h5ad
notebooks/20_scgpt_bcl6.ipynb        -> data/processed/perturbations/scgpt_bcl6_native.h5ad
notebooks/30_geneformer_bcl6.ipynb   -> data/processed/perturbations/geneformer_bcl6_native.npz
```

Run order for one model (CellOracle shown):

```bash
conda activate bcl6-celloracle
jupyter notebook notebooks/10_celloracle_bcl6.ipynb     # produces native h5ad
conda activate bcl6-recovery
python scripts/07_run_celloracle_bcl6.py \
    --config configs/config.yaml                         # native -> standard CSV
python scripts/05_compare_perturbation_models.py \
    --config configs/config.yaml                         # score and compare
```

To enable a runner inside the orchestrator (`scripts/06_run_model_perturbations.py`
and `scripts/10_run_end_to_end_pipeline.py`), flip `enabled: true` for that
runner in `configs/config.yaml`.

### Bring-your-own-result

If you ran any of these models elsewhere and already have a native artefact,
save it in one of these forms (see `scripts/converters/*` for the auto-detect
logic):

| Model | Acceptable native artefact |
|-------|----------------------------|
| CellOracle | `.h5ad` with layers `imputed_count` & `simulated_count`, or `.npz` with `original` / `perturbed` / `gene_names` / `cell_ids`, or pickled Oracle (`.pkl`) post `simulate_shift` |
| scGPT | `.h5ad` with `.layers["perturbed"]` (or `.X` = perturbed), or `.npz` with `perturbed` / `gene_names` / `cell_ids`, or genes x cells `.csv` |
| Geneformer | `.npz` with `cell_ids` / `original_embeddings` / `perturbed_embeddings` |

Then point `model_runners.<model>.native_output_path` at the file (or pass
`--native-output PATH` on the CLI) and run the matching `scripts/0X_run_*.py`.

### Notes on cross-model comparability

- CellOracle outputs a delta in the same scale as the input matrix
  (`delta_expression` format).
- scGPT outputs predicted post-perturbation expression (`perturbed_expression`).
- Geneformer outputs only embedding shifts. We project them back to expression
  space via ridge regression fitted on the exported expression matrix
  (`scripts/converters/geneformer_native.py`). Tune `ridge_alpha` in
  `model_runners.geneformer` if the projection looks unstable.
- All three feed the same `delta_recovery_probability` /
  `recovery_shift_score` evaluator, so a per-model bias in absolute scale is
  partially absorbed by the per-cell normalisation built into those metrics —
  but absolute deltas should still be interpreted with care.

You can also run just the model collection step (which calls the converters):

```bash
python scripts/06_run_model_perturbations.py --config configs/config.yaml
```

## Important caution

The ChIP file provides BCL6 target genes but not the regulatory direction. The first-pass perturbation script infers target direction using correlation with BCL6 expression. This is useful for exploratory analysis, but a stronger mechanistic version should use CellOracle/GRN-based perturbation.
