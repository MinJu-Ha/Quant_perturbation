"""Shared helpers for converting model-specific perturbation outputs into the
project's standard genes x cells CSV format.

The standard CSV is what scripts/05_compare_perturbation_models.py consumes:

    gene,cell_1,cell_2,...
    GeneA,0.12,-0.03,...
    GeneB,-0.31,0.22,...

Rows are genes. The first column is the gene symbol. Columns are cell IDs that
match data/processed/Cardiomyocyte_meta.csv.
"""
from pathlib import Path

import numpy as np
import pandas as pd


def read_expression_matrix(path):
    expr = pd.read_csv(path, index_col=0)
    expr.index = expr.index.astype(str)
    expr.columns = expr.columns.astype(str)
    return expr


def write_genes_by_cells_csv(matrix, path, index_label="gene"):
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    matrix = matrix.copy()
    matrix.index = matrix.index.astype(str)
    matrix.columns = matrix.columns.astype(str)
    matrix.index.name = index_label
    matrix.to_csv(path)


def to_genes_by_cells(matrix_df, gene_axis):
    if gene_axis == "rows":
        return matrix_df
    if gene_axis == "columns":
        return matrix_df.T
    raise ValueError(f"gene_axis must be 'rows' or 'columns', got {gene_axis!r}")


def align_to_reference_csv(matrix_df, reference_csv_path, fill_value=0.0):
    """Reindex a genes x cells DataFrame to match a reference expression matrix's
    gene rows and cell columns. Missing cells/genes are filled with fill_value.

    Use fill_value=0.0 with delta_expression outputs, so unmodelled cells get
    a no-effect delta when the project evaluator merges them in.
    """
    reference = read_expression_matrix(reference_csv_path)
    aligned = matrix_df.reindex(
        index=reference.index,
        columns=reference.columns,
        fill_value=fill_value,
    )
    return aligned


def safe_to_dense(arr):
    if hasattr(arr, "toarray"):
        return np.asarray(arr.toarray())
    return np.asarray(arr)
