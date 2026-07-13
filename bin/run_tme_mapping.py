#!/usr/bin/env python
"""
TME mapping: neighborhood enrichment analysis on deconvolved cell type proportions.
Determines which cell type pairs co-localize spatially more than expected by chance.

Input: deconvolution output .h5ad with cell type proportions in .obsm
Output: .h5ad with neighborhood enrichment scores + summary text file
"""
import argparse
import numpy as np
import pandas as pd
import scanpy as sc
import squidpy as sq

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--h5ad", required=True, help="Path to deconvolution .h5ad")
    parser.add_argument("--celltype_col", default="celltype_major")
    parser.add_argument("--n_perms", type=int, default=1000)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    adata = sc.read_h5ad(args.h5ad)

    # assign each spot a dominant cell type based on highest estimated abundance,
    # so we can run categorical neighborhood enrichment (squidpy's standard approach)
    props = adata.obsm["q05_cell_abundance_w_sf"]
    props.columns = [c.replace("q05cell_abundance_w_sf_", "") for c in props.columns]
    adata.obs["dominant_celltype"] = props.idxmax(axis=1).values
    adata.obs["dominant_celltype"] = adata.obs["dominant_celltype"].astype("category")

    print("Dominant cell type distribution across spots:")
    print(adata.obs["dominant_celltype"].value_counts())

    sq.gr.spatial_neighbors(adata)
    sq.gr.nhood_enrichment(adata, cluster_key="dominant_celltype", n_perms=args.n_perms, seed=0)

    zscore = adata.uns["dominant_celltype_nhood_enrichment"]["zscore"]
    categories = adata.obs["dominant_celltype"].cat.categories
    enrichment_df = pd.DataFrame(zscore, index=categories, columns=categories)

    adata.write_h5ad(args.out)

    with open(args.out.replace(".h5ad", "_tme_summary.txt"), "w") as f:
        f.write("Dominant cell type distribution:\n")
        f.write(adata.obs["dominant_celltype"].value_counts().to_string())
        f.write("\n\nNeighborhood enrichment z-scores (positive = co-localize more than chance):\n")
        f.write(enrichment_df.round(2).to_string())
        f.write("\n")

    top_pair_idx = np.unravel_index(np.argmax(zscore - np.diag(np.diag(zscore))), zscore.shape)
    print(f"TME mapping done: strongest co-localization = {categories[top_pair_idx[0]]} <-> {categories[top_pair_idx[1]]}, z={zscore[top_pair_idx]:.2f}")

if __name__ == "__main__":
    main()
