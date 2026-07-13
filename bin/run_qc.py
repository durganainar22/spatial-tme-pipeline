#!/usr/bin/env python
"""
Step 1 of pipeline: Quality control on Visium spatial data.
Input: AnnData .h5ad file (spots x genes, with spatial coords in .obsm)
Output: filtered AnnData object (.h5ad) + a QC summary text file
"""
import argparse
import scanpy as sc

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--h5ad", required=True, help="Path to input .h5ad file")
    parser.add_argument("--min_genes", type=int, default=200)
    parser.add_argument("--max_mito_pct", type=float, default=20.0)
    parser.add_argument("--out", required=True, help="Output .h5ad path")
    args = parser.parse_args()

    adata = sc.read_h5ad(args.h5ad)
    adata.var_names_make_unique()

    n_spots_before = adata.n_obs

    adata.var["mt"] = adata.var_names.str.startswith("MT-")
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True, percent_top=None)

    sc.pp.filter_cells(adata, min_genes=args.min_genes)
    adata = adata[adata.obs["pct_counts_mt"] < args.max_mito_pct].copy()

    n_spots_after = adata.n_obs

    adata.write_h5ad(args.out)

    with open(args.out.replace(".h5ad", "_qc_summary.txt"), "w") as f:
        f.write(f"Spots before QC: {n_spots_before}\n")
        f.write(f"Spots after QC: {n_spots_after}\n")
        f.write(f"Filters: min_genes={args.min_genes}, max_mito_pct={args.max_mito_pct}\n")

    print(f"QC done: {n_spots_before} -> {n_spots_after} spots")

if __name__ == "__main__":
    main()
