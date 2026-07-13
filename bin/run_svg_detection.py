#!/usr/bin/env python
"""
Step 3 of pipeline: Detect spatially variable genes using Moran's I statistic.
Input: normalized AnnData (.h5ad) with spatial coords in .obsm['spatial']
Output: AnnData with Moran's I scores added + a summary text file
"""
import argparse
import scanpy as sc
import squidpy as sq

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--h5ad", required=True, help="Path to normalized .h5ad file")
    parser.add_argument("--n_top_svgs", type=int, default=50)
    parser.add_argument("--out", required=True, help="Output .h5ad path")
    args = parser.parse_args()

    adata = sc.read_h5ad(args.h5ad)

    # build spatial neighbor graph -- required before Moran's I
    sq.gr.spatial_neighbors(adata)

    # restrict to HVGs to keep this fast and biologically focused
    sq.gr.spatial_autocorr(adata, mode="moran", genes=adata.var_names[adata.var["highly_variable"]], n_perms=100, n_jobs=1)

    moran_results = adata.uns["moranI"]
    top_svgs = moran_results.head(args.n_top_svgs)

    adata.write_h5ad(args.out)

    with open(args.out.replace(".h5ad", "_svg_summary.txt"), "w") as f:
        f.write(f"Genes tested (HVGs): {int(adata.var['highly_variable'].sum())}\n")
        f.write(f"Top {args.n_top_svgs} spatially variable genes:\n")
        f.write(top_svgs.to_string())
        f.write("\n")

    print(f"SVG detection done: top gene = {top_svgs.index[0]}, Moran's I = {top_svgs['I'].iloc[0]:.4f}")

if __name__ == "__main__":
    main()
