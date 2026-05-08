# BCL6 Recovery-State Perturbation Project

Quantitatively test whether **BCL6 in silico perturbation** shifts post-LVAD
cardiomyocytes toward a **recovery-like transcriptional state**, and compare
the answer across **CellOracle**, **scGPT**, and **Geneformer**.

This repo only owns the **shared scoring/quantification path**. Each model is
run upstream (in its own pipeline / notebook). Drop the result into
`models/<Name>/results/` and the scripts below pick it up.

## Inputs and exports

```text
data/raw/post_cardiomyocyte.rds
data/raw/Bcl6_ChIP_target_P0.05.csv
```

Step 1 exports them into:

```text
data/processed/Cardiomyocyte_meta.csv
data/processed/Cardiomyocyte_SCT_scaled.txt   # SCT scale.data, BCL6 + ChIP genes
data/processed/Cardiomyocyte_RNA_counts.txt   # RNA counts slot, BCL6 + ChIP genes
```

The RNA counts file is what scGPT / Geneformer should consume. SCT_scaled is
what the recovery classifier and centroid metrics use.

## Labels

```text
recovered:     Rpre, Rpost
non-recovered: NRpre, NRpost
reference:     Donor
```

Edit them in `configs/config.yaml -> metadata`.

## Per-model output layout

Each model lives under its own folder. Drop either a native artefact or a
ready-made standard CSV in `results/`:

```text
models/CellOracle/results/
    standard_perturbation.csv             # gene x cell, delta_expression
    celloracle_bcl6_native.h5ad           # OR: AnnData with imputed_count + simulated_count layers
models/scGPT/results/
    standard_perturbation.csv             # gene x cell, perturbed_expression
    scgpt_bcl6_native.h5ad                # OR: AnnData with .layers["perturbed"]
models/Geneformer/results/
    geneformer_bcl6_native.npz            # cell_ids + original_embeddings + perturbed_embeddings
                                          # script 03 projects to expression via ridge regression
```

Add / remove models or change paths in `configs/config.yaml -> model_results.models[]`.

## Scripts (run in order)

```text
scripts/01_export_from_seurat.R                      # rds -> meta + SCT_scaled + RNA counts
scripts/02_train_recovery_classifier.py              # logistic recovered-vs-nonrecovered classifier
scripts/03_check_and_transform_model_results.py      # native -> standard CSV + manifest
scripts/04_quantify_model_results.py                 # apply classifier + centroid metrics per model
scripts/05_compare_and_plot_model_results.py         # pairwise concordance + box plots
```

Helpers:

```text
scripts/replay_celloracle.py    # one-off: replay simulate_shift on an existing oracle
scripts/converters/             # native -> standard converters used by script 03
```

## Quick start

```bash
# 1. Data prep (R, slow on 5 GB rds)
Rscript scripts/01_export_from_seurat.R

# 2. Train recovery classifier
python scripts/02_train_recovery_classifier.py --config configs/config.yaml

# 3. Drop each model's output into models/<Name>/results/, then:
python scripts/03_check_and_transform_model_results.py --config configs/config.yaml
python scripts/04_quantify_model_results.py           --config configs/config.yaml
python scripts/05_compare_and_plot_model_results.py   --config configs/config.yaml
```

Missing models are skipped automatically by step 3 (see
`results/tables/model_result_manifest.csv`).

## CellOracle replay (shortcut)

If you already ran the BCL6 KO simulate_shift in a CellOracle notebook but
did not save the post-perturbation oracle, replay just that step:

```bash
conda activate celloracle
python scripts/replay_celloracle.py \
    --oracle <path>/Cardiomyocyte2.celloracle.oracle \
    --links  <path>/links_cardiomyocyte2_cm.subtype.celloracle.links \
    --target BCL6 \
    --out    models/CellOracle/results/celloracle_bcl6_native.h5ad
```

## Quantification metrics

Per cell:

```text
delta_recovery_probability = perturbed_recovery_probability - original_recovery_probability
recovery_shift_score       = dist(orig, recovered_centroid) - dist(perturbed, recovered_centroid)
nonrecovery_escape_score   = dist(perturbed, nonrec_centroid) - dist(orig, nonrec_centroid)
```

`> 0` means the perturbation moves the cell toward the recovered state.
Outputs: `results/tables/model_quantification_*.csv`,
`results/figures/model_quantification_*.png`.

## Caution

The ChIP file lists BCL6 target genes but not regulatory direction. The
recovery classifier learns a data-driven sign per target gene; the
perturbation models themselves provide the direction. The quantification is
only as meaningful as each model's perturbation prediction.
