#!/usr/bin/env python
"""
Cell2location Stage 2: Map reference cell type signatures onto spatial data
to estimate cell type abundance per Visium spot.

Input: spatial .h5ad (raw counts) + reference_signatures.h5ad (from Stage 1)
Output: spatial .h5ad with per-spot cell type proportions in .obsm
"""
import argparse
import numpy as np
import scanpy as sc
import anndata as ad
import cell2location
from cell2location.models import Cell2location

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--spatial_h5ad", required=True)
    parser.add_argument("--ref_signatures", required=True)
    parser.add_argument("--celltype_col", default="celltype_major")
    parser.add_argument("--max_epochs", type=int, default=500)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    adata_ref = sc.read_h5ad(args.ref_signatures)

    # reference uses Ensembl IDs as var_names; spatial data uses gene symbols.
    # remap reference to symbols via CELLxGENE's feature_name column, dropping
    # genes with no real symbol (where feature_name == Ensembl ID) and deduping.
    adata_ref.var["gene_symbol"] = adata_ref.var["feature_name"].astype(str)
    has_real_symbol = adata_ref.var["gene_symbol"] != adata_ref.var_names.astype(str)
    adata_ref = adata_ref[:, has_real_symbol].copy()
    adata_ref.var_names = adata_ref.var["gene_symbol"]
    adata_ref = adata_ref[:, ~adata_ref.var_names.duplicated()].copy()
    print(f"Reference after symbol remap + dedup: {adata_ref.n_vars} genes")

    # extract per-cell-type mean expression signatures from the trained reference model
    inf_aver = adata_ref.varm["means_per_cluster_mu_fg"].copy()
    inf_aver.index = adata_ref.var_names
    inf_aver.columns = adata_ref.obs[args.celltype_col].cat.categories

    adata_vis = sc.read_h5ad(args.spatial_h5ad)
    print(f"Loaded spatial data: {adata_vis.n_obs} spots x {adata_vis.n_vars} genes")

    if adata_vis.raw is not None:
        adata_vis = ad.AnnData(X=adata_vis.raw.X, obs=adata_vis.obs, var=adata_vis.raw.var, obsm=adata_vis.obsm)
    elif "counts" in adata_vis.layers:
        adata_vis.X = adata_vis.layers["counts"].copy()

    shared_genes = [g for g in inf_aver.index if g in adata_vis.var_names]
    adata_vis = adata_vis[:, shared_genes].copy()
    inf_aver = inf_aver.loc[shared_genes]
    print(f"Shared genes for mapping: {len(shared_genes)}")

    if len(shared_genes) < 100:
        raise ValueError(f"Only {len(shared_genes)} shared genes found -- something is likely still wrong with gene ID matching")

    Cell2location.setup_anndata(adata=adata_vis)

    mod = Cell2location(
        adata_vis,
        cell_state_df=inf_aver,
        N_cells_per_location=30,
        detection_alpha=20,
    )

    mod.train(max_epochs=args.max_epochs, accelerator="gpu", train_size=1)

    adata_vis = mod.export_posterior(
        adata_vis, sample_kwargs={"num_samples": 1000, "batch_size": adata_vis.n_obs, "accelerator": "gpu"}
    )

    adata_vis.write_h5ad(args.out)
    print(f"Deconvolution done: cell type proportions added for {adata_vis.n_obs} spots")

if __name__ == "__main__":
    main()
