"""Convert CellOracle perturbation outputs into the project's standard
genes x cells delta-expression CSV.

CellOracle produces a perturbed expression matrix as part of its
`oracle.simulate_shift()` step. The native artefact this converter accepts is
auto-detected by extension:

  *.h5ad : AnnData with layers
            - "imputed_count"   (cells x genes) original imputed expression
            - "simulated_count" (cells x genes) post-perturbation expression
  *.npz  : numpy archive with arrays
            - "original"   (cells x genes)
            - "perturbed"  (cells x genes)
            - "gene_names" (genes,)
            - "cell_ids"   (cells,)
  *.pkl  : pickled CellOracle Oracle object that has already had simulate_shift
            called on it. The same `imputed_count` / `simulated_count` layers
            are read from `oracle.adata`.

Anyone who already has a CellOracle Oracle from their own pipeline can just
pickle it and feed it here — they do not need our notebook.
"""
from pathlib import Path
import pickle

import numpy as np
import pandas as pd

from .common import write_genes_by_cells_csv


def _layers_to_delta(adata):
    if "imputed_count" not in adata.layers or "simulated_count" not in adata.layers:
        raise ValueError(
            "AnnData/Oracle must contain layers 'imputed_count' and "
            "'simulated_count' produced by CellOracle.simulate_shift()."
        )
    imputed = np.asarray(adata.layers["imputed_count"])
    simulated = np.asarray(adata.layers["simulated_count"])
    delta = simulated - imputed  # cells x genes
    return pd.DataFrame(
        delta.T,
        index=np.asarray(adata.var_names).astype(str),
        columns=np.asarray(adata.obs_names).astype(str),
    )


def _from_anndata(path):
    import anndata as ad
    adata = ad.read_h5ad(path)
    return _layers_to_delta(adata)


def _from_npz(path):
    z = np.load(path, allow_pickle=True)
    delta = np.asarray(z["perturbed"]) - np.asarray(z["original"])
    gene_names = np.asarray(z["gene_names"]).astype(str)
    cell_ids = np.asarray(z["cell_ids"]).astype(str)
    return pd.DataFrame(delta.T, index=gene_names, columns=cell_ids)


def _from_pickle(path):
    with open(path, "rb") as f:
        oracle = pickle.load(f)
    return _layers_to_delta(oracle.adata)


def convert(native_path, output_path):
    native_path = Path(native_path)
    suffix = native_path.suffix.lower()
    if suffix == ".h5ad":
        delta_df = _from_anndata(native_path)
    elif suffix == ".npz":
        delta_df = _from_npz(native_path)
    elif suffix in {".pkl", ".pickle", ".oracle"}:
        delta_df = _from_pickle(native_path)
    else:
        raise ValueError(
            f"Unsupported CellOracle native_path extension: {suffix}. "
            "Expected .h5ad, .npz, .pkl, .pickle, or .oracle."
        )
    write_genes_by_cells_csv(delta_df, output_path)
    return delta_df.shape
