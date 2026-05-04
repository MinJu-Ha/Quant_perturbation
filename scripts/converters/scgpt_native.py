"""Convert scGPT perturbation outputs into the project's standard
genes x cells perturbed-expression CSV.

scGPT's perturbation utilities (zero-shot KO via input masking, or fine-tuned
GEARS-style prediction) all eventually produce a cells x genes matrix of
predicted expression. This converter normalises whatever artefact the user
saved into the project's standard format.

Auto-detected native input formats:

  *.h5ad : AnnData. Looked up in this order:
              - .layers["perturbed"]   (preferred; cells x genes)
              - .layers["X_perturbed"] (alias)
              - .X                     (assumed perturbed; cells x genes)
  *.npz  : numpy archive with arrays
              - "perturbed"  (cells x genes)   required
              - "gene_names" (genes,)          required
              - "cell_ids"   (cells,)          required
              - "original"   (cells x genes)   optional
  *.csv  : already a genes x cells perturbed-expression matrix; this converter
            just copies/renames it. Useful when the user has a raw export.

Anyone who already ran scGPT independently can pickle/save into any of these
shapes — they do not need our notebook.
"""
from pathlib import Path

import numpy as np
import pandas as pd

from .common import write_genes_by_cells_csv, safe_to_dense


def _from_anndata(path):
    import anndata as ad
    adata = ad.read_h5ad(path)
    layer_keys = ["perturbed", "X_perturbed"]
    perturbed = None
    for key in layer_keys:
        if key in adata.layers:
            perturbed = safe_to_dense(adata.layers[key])
            break
    if perturbed is None:
        perturbed = safe_to_dense(adata.X)

    df = pd.DataFrame(
        perturbed.T,
        index=np.asarray(adata.var_names).astype(str),
        columns=np.asarray(adata.obs_names).astype(str),
    )
    return df


def _from_npz(path):
    z = np.load(path, allow_pickle=True)
    perturbed = np.asarray(z["perturbed"])
    gene_names = np.asarray(z["gene_names"]).astype(str)
    cell_ids = np.asarray(z["cell_ids"]).astype(str)
    return pd.DataFrame(perturbed.T, index=gene_names, columns=cell_ids)


def _from_csv(path):
    df = pd.read_csv(path, index_col=0)
    df.index = df.index.astype(str)
    df.columns = df.columns.astype(str)
    return df


def convert(native_path, output_path):
    native_path = Path(native_path)
    suffix = native_path.suffix.lower()
    if suffix == ".h5ad":
        df = _from_anndata(native_path)
    elif suffix == ".npz":
        df = _from_npz(native_path)
    elif suffix in {".csv", ".tsv", ".txt"}:
        df = _from_csv(native_path)
    else:
        raise ValueError(
            f"Unsupported scGPT native_path extension: {suffix}. "
            "Expected .h5ad, .npz, or .csv."
        )
    write_genes_by_cells_csv(df, output_path)
    return df.shape
