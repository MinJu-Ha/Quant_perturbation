"""Convert Geneformer InSilicoPerturber outputs into the project's standard
genes x cells perturbed-expression CSV.

Geneformer is fundamentally an embedding model: its native perturbation
artefact is a cell-state shift in embedding space (n_cells x emb_dim), not a
gene-level expression change. To plug into our shared evaluator (which expects
genes x cells), we project embedding deltas back into expression space using a
ridge-regression linear map fitted on the same dataset's original embeddings
and original expression.

This is a deliberate, principled approximation:

    expr_orig  ≈  embedding_orig  @ W + b      (fit by ridge)
    expr_pert  =  embedding_pert  @ W + b      (apply the same map)

Auto-detected native input formats:

  *.npz  : numpy archive with arrays
            - "cell_ids"            (cells,)         required
            - "original_embeddings" (cells x emb_dim) required
            - "perturbed_embeddings"(cells x emb_dim) required

The original expression matrix (the same `Cardiomyocyte_SCT_scaled.txt` exported
by scripts/01_export_from_seurat.R) is required for the projection fit. Pass it
via the `original_expression_csv` argument or rely on the default in
`scripts/09_run_geneformer_bcl6.py`.

If a user already has Geneformer embeddings from their own pipeline, saving
them as a single .npz with the three arrays above is enough to plug into this
converter.
"""
from pathlib import Path

import numpy as np
import pandas as pd

from .common import read_expression_matrix, write_genes_by_cells_csv


def _load_embeddings(native_path):
    native_path = Path(native_path)
    suffix = native_path.suffix.lower()
    if suffix == ".npz":
        z = np.load(native_path, allow_pickle=True)
        cell_ids = np.asarray(z["cell_ids"]).astype(str)
        orig = np.asarray(z["original_embeddings"], dtype=float)
        pert = np.asarray(z["perturbed_embeddings"], dtype=float)
        return cell_ids, orig, pert
    raise ValueError(
        f"Unsupported Geneformer native_path extension: {suffix}. Expected .npz."
    )


def _fit_ridge_projection(emb, expr_cells_by_genes, ridge_alpha=1.0):
    """Solve expr ≈ emb @ W + b via ridge regression.

    emb: (n_cells, emb_dim)
    expr_cells_by_genes: (n_cells, n_genes)
    Returns (W, b) with shapes (emb_dim, n_genes) and (n_genes,).
    """
    emb_mean = emb.mean(axis=0)
    expr_mean = expr_cells_by_genes.mean(axis=0)
    Ec = emb - emb_mean
    Yc = expr_cells_by_genes - expr_mean
    n_emb = Ec.shape[1]
    A = Ec.T @ Ec + ridge_alpha * np.eye(n_emb)
    B = Ec.T @ Yc
    W = np.linalg.solve(A, B)
    b = expr_mean - emb_mean @ W
    return W, b


def convert(native_path, output_path, original_expression_csv, ridge_alpha=1.0):
    cell_ids, orig_emb, pert_emb = _load_embeddings(native_path)
    expr = read_expression_matrix(original_expression_csv)

    common_cells = pd.Index(cell_ids).intersection(expr.columns)
    if len(common_cells) == 0:
        raise ValueError(
            "No overlap between Geneformer cell_ids and the exported expression "
            "matrix columns. Ensure cell_ids in the .npz match "
            "data/processed/Cardiomyocyte_SCT_scaled.txt column names."
        )

    cell_to_idx = {cid: i for i, cid in enumerate(cell_ids)}
    keep = np.array([cell_to_idx[c] for c in common_cells])
    orig_emb_aligned = orig_emb[keep]
    pert_emb_aligned = pert_emb[keep]
    expr_aligned = expr.loc[:, common_cells].values.T  # cells x genes

    W, b = _fit_ridge_projection(orig_emb_aligned, expr_aligned, ridge_alpha=ridge_alpha)

    pert_expr = pert_emb_aligned @ W + b  # cells x genes
    df = pd.DataFrame(
        pert_expr.T,
        index=expr.index.astype(str),
        columns=common_cells.astype(str),
    )
    write_genes_by_cells_csv(df, output_path)
    return df.shape
