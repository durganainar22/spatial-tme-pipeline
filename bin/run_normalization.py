#!/usr/bin/env python
"""
Step 2 of pipeline: Normalize counts and select highly variable genes.
Input: QC-filtered AnnData (.h5ad)
Output: normalized AnnData (.h5ad) with HVGs flagged + a summary text file
"""
import argparse
import scanpy as sc

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--h5ad", required=True, help="Path to QC-filtered .h5ad file")
    parser.add_argument("--n_top_genes", type=int, default=2000)
    parser.add_argument("--out", required=True, help="Output .h5ad path")
    args = parser.parse_args()

    adata = sc.read_h5ad(args.h5ad)

    # keep raw counts before normalizing, common practice for downstream methods
    adata.layers["counts"] = adata.X.copy()

    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

    sc.pp.highly_variable_genes(adata, n_top_genes=args.n_top_genes, flavor="seurat")

    n_hvg = int(adata.var["highly_variable"].sum())

    adata.write_h5ad(args.out)

    with open(args.out.replace(".h5ad", "_norm_summary.txt"), "w") as f:
        f.write(f"Total genes: {adata.n_vars}\n")
        f.write(f"Highly variable genes selected: {n_hvg}\n")
        f.write(f"n_top_genes requested: {args.n_top_genes}\n")

    print(f"Normalization done: {n_hvg} HVGs selected out of {adata.n_vars} genes")

if __name__ == "__main__":
    main()
