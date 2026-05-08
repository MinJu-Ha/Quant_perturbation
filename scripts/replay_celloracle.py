"""Reproduce the BCL6 KO simulate_shift step on a memory-light path.

Loading the user's saved 12 GB pickled Oracle OOM-kills on a 24 GB Mac. This
script instead rebuilds the Oracle from the much smaller `cardiomyocyte2.h5ad`
plus the existing `.links` file, mirroring the preprocessing in
`Minju_bcl6/1_Network_analysis_cardio2.ipynb`:

  * relabel `obs['condition']` (Donor/preLVAD HF/NRpost/Rpost),
  * subsample to 20000 cells (random_state=123),
  * subset to vst.variable genes,
  * Oracle.import_anndata_as_raw_count(cluster='cm.subtype', embedding='X_umap'),
  * perform_PCA, knn_imputation with k = round(0.025 * n_cells),
  * load existing links, fit_GRN_for_simulation, simulate_shift({BCL6: 0.0}),
  * write a compact .h5ad with imputed_count and simulated_count layers.

The output is what `scripts/converters/celloracle_native.py` consumes via
`scripts/03_check_and_transform_model_results.py`.
"""
from pathlib import Path
import argparse
import sys
import types

import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad


def _install_pandas_legacy_shim():
    """Re-expose `pandas.core.indexes.numeric` for unpickling artefacts saved
    under pandas <2.0. The module was merged into base in 2.x; CellOracle's
    saved Links pickle still references the old path.
    """
    try:
        import pandas.core.indexes.numeric  # noqa: F401
        return
    except ImportError:
        pass
    from pandas import Index, Int64Dtype, Float64Dtype  # noqa: F401
    shim = types.ModuleType("pandas.core.indexes.numeric")
    shim.NumericIndex = Index
    shim.Int64Index = Index
    shim.UInt64Index = Index
    shim.Float64Index = Index
    sys.modules["pandas.core.indexes.numeric"] = shim


_install_pandas_legacy_shim()
import celloracle as co


CONDITION_RELABEL = {
    "Donor": "Donor",
    "NRpre": "preLVAD HF",
    "Rpre": "preLVAD HF",
    "NRpost": "NRpost",
    "Rpost": "Rpost",
}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--h5ad", default="data/processed/celloracle_input/Minju_bcl6/cardiomyocyte2.h5ad")
    ap.add_argument("--links", default="/Users/minjuha/Downloads/Minju_bcl6/links_cardiomyocyte2_cm.subtype.celloracle.links")
    ap.add_argument("--cluster-col", default="cm.subtype")
    ap.add_argument("--embedding", default="X_umap")
    ap.add_argument("--target", default="BCL6")
    ap.add_argument("--n-propagation", type=int, default=3)
    ap.add_argument("--alpha", type=float, default=10.0)
    ap.add_argument("--n-cells-downsample", type=int, default=20000)
    ap.add_argument("--downsample-seed", type=int, default=123)
    ap.add_argument("--n-jobs", type=int, default=4)
    ap.add_argument("--out", default="models/CellOracle/results/celloracle_bcl6_native.h5ad")
    args = ap.parse_args()

    h5ad_path = Path(args.h5ad)
    links_path = Path(args.links)
    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    print(f"[1/8] celloracle {co.__version__}")

    print(f"[2/8] Reading h5ad: {h5ad_path}")
    adata = sc.read_h5ad(h5ad_path)
    print(f"      shape {adata.shape}, layers {list(adata.layers.keys())}, obsm {list(adata.obsm.keys())}")

    if "raw_count" in adata.layers:
        adata.X = adata.layers["raw_count"].copy()
    print(f"      var columns {list(adata.var.columns)[:8]}")
    if "vst.variable" in adata.var.columns:
        adata.var["variable_gene"] = adata.var["vst.variable"].astype(bool)
    elif "variable_gene" not in adata.var.columns:
        raise SystemExit("h5ad var has neither 'vst.variable' nor 'variable_gene'.")

    print(f"[3/8] Relabel condition (per notebook 1)")
    adata.obs["condition"] = adata.obs["condition"].astype(str).replace(CONDITION_RELABEL)
    print(f"      condition counts: {adata.obs['condition'].value_counts().to_dict()}")

    if adata.shape[0] > args.n_cells_downsample:
        print(f"[4/8] Subsample to {args.n_cells_downsample} cells (seed={args.downsample_seed})")
        sc.pp.subsample(adata, n_obs=args.n_cells_downsample, random_state=args.downsample_seed)
    print(f"      cells now {adata.shape[0]}")

    print(f"[5/8] Subset to variable genes")
    var_genes = adata.var.index[adata.var["variable_gene"]].values
    adata = adata[:, var_genes].copy()
    print(f"      shape after var-gene filter {adata.shape}")
    if args.target not in adata.var_names:
        raise SystemExit(f"Target gene {args.target} not in variable-gene set.")

    print(f"[6/8] Build Oracle + base GRN + PCA + KNN imputation")
    oracle = co.Oracle()
    oracle.import_anndata_as_raw_count(
        adata=adata,
        cluster_column_name=args.cluster_col,
        embedding_name=args.embedding,
    )
    base_GRN = co.data.load_human_promoter_base_GRN()
    oracle.import_TF_data(TF_info_matrix=base_GRN)
    oracle.perform_PCA()
    cum = np.cumsum(oracle.pca.explained_variance_ratio_)
    diff_diff = np.diff(np.diff(cum) > 0.002)
    knee = np.where(diff_diff)[0]
    n_comps = int(knee[0]) if knee.size else 50
    n_comps = min(n_comps, 50)
    n_cells = oracle.adata.shape[0]
    k = int(0.025 * n_cells)
    print(f"      n_comps={n_comps}, k={k}")
    oracle.knn_imputation(
        n_pca_dims=n_comps,
        k=k,
        balanced=True,
        b_sight=k * 8,
        b_maxl=k * 4,
        n_jobs=args.n_jobs,
    )

    print(f"[7/8] Load links + fit GRN for simulation")
    links = co.load_hdf5(str(links_path))
    # The saved Links object already has filtered_links computed from the original
    # network analysis (threshold_number=2000). Re-running filter_links() with the
    # current celloracle's defaults rebuilds it differently and drops BCL6 from
    # most clusters' TF dict — so we keep the saved filter as-is.
    oracle.get_cluster_specific_TFdict_from_Links(links_object=links)
    oracle.fit_GRN_for_simulation(alpha=args.alpha, use_cluster_specific_TFdict=True)

    print(f"[8/8] simulate_shift {args.target} -> 0.0 (n_propagation={args.n_propagation})")
    oracle.simulate_shift(perturb_condition={args.target: 0.0}, n_propagation=args.n_propagation)

    layers = oracle.adata.layers
    if "imputed_count" not in layers or "simulated_count" not in layers:
        raise SystemExit(f"Missing expected layers, have {list(layers.keys())}")

    keep = ad.AnnData(
        X=np.asarray(layers["imputed_count"]).copy(),
        obs=oracle.adata.obs[[args.cluster_col, "condition"]].copy(),
        var=pd.DataFrame(index=oracle.adata.var_names.astype(str)),
    )
    keep.layers["imputed_count"] = np.asarray(layers["imputed_count"]).copy()
    keep.layers["simulated_count"] = np.asarray(layers["simulated_count"]).copy()
    keep.uns["perturbation"] = {
        "target_gene": args.target,
        "kind": "celloracle_simulate_shift_fresh_build",
        "n_propagation": int(args.n_propagation),
        "alpha": float(args.alpha),
        "n_cells_downsample": int(args.n_cells_downsample),
        "downsample_seed": int(args.downsample_seed),
        "n_pca_dims": int(n_comps),
        "k": int(k),
    }
    keep.write_h5ad(out_path, compression="gzip")
    print(f"Done. Wrote {out_path} ({out_path.stat().st_size/1e6:.1f} MB)")


if __name__ == "__main__":
    main()
