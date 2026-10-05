#!/usr/bin/env python
"""
Cell2location Stage 1: Estimate reference cell type expression signatures
from single-cell reference data (Wu et al. 2021 breast cancer atlas).
Input: single-cell reference .h5ad with a cell type annotation column
Output: reference signatures (.h5ad model output) used later for spatial mapping
"""
import argparse
import scanpy as sc
import anndata as ad
import cell2location
from cell2location.utils.filtering import filter_genes
from cell2location.models import RegressionModel

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--ref_h5ad", required=True, help="Path to single-cell reference .h5ad")
    parser.add_argument("--celltype_col", default="celltype_major", help="obs column with cell type labels")
    parser.add_argument("--batch_col", default="donor_id", help="obs column for batch/donor")
    parser.add_argument("--max_epochs", type=int, default=250)
    parser.add_argument("--accelerator", default="gpu", help="gpu, cpu or auto")
    parser.add_argument("--out", required=True, help="Output directory for reference model + signatures")
    args = parser.parse_args()

    adata_ref = sc.read_h5ad(args.ref_h5ad)
    print(f"Loaded reference: {adata_ref.n_obs} cells x {adata_ref.n_vars} genes")

    # cell2location requires raw integer counts, not normalized/log data.
    # This dataset stores raw counts in .raw -- pull them into .X before anything else.
    if adata_ref.raw is not None:
        adata_ref = ad.AnnData(
            X=adata_ref.raw.X,
            obs=adata_ref.obs,
            var=adata_ref.raw.var,
        )
        print(f"Restored raw counts from .raw: X max = {adata_ref.X.max()}")
    else:
        raise ValueError("No .raw counts found -- cannot proceed without raw integer counts")

    # standard cell2location gene filtering -- removes lowly expressed / uninformative genes
    selected = filter_genes(adata_ref, cell_count_cutoff=5, cell_percentage_cutoff2=0.03, nonz_mean_cutoff=1.12)
    adata_ref = adata_ref[:, selected].copy()
    print(f"After filtering: {adata_ref.n_vars} genes retained")

    RegressionModel.setup_anndata(
        adata=adata_ref,
        batch_key=args.batch_col,
        labels_key=args.celltype_col,
    )

    mod = RegressionModel(adata_ref)
    mod.train(max_epochs=args.max_epochs, accelerator=args.accelerator)

    adata_ref = mod.export_posterior(
        adata_ref, sample_kwargs={"num_samples": 1000, "batch_size": 2500, "accelerator": args.accelerator}
    )

    mod.save(f"{args.out}/reference_model", overwrite=True)
    adata_ref.write_h5ad(f"{args.out}/reference_signatures.h5ad")

    n_celltypes = adata_ref.obs[args.celltype_col].nunique()
    print(f"Reference signature estimation done: {n_celltypes} cell types, {adata_ref.n_vars} genes")

if __name__ == "__main__":
    main()
