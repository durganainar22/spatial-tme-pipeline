#!/usr/bin/env python
"""
Final pipeline stage: visualize spatial cell type distribution and
neighborhood enrichment results.

Input: tme.h5ad (from TME_MAPPING stage)
Output: spatial scatter plot, enrichment heatmap, summary report
"""
import argparse
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--h5ad", required=True)
    parser.add_argument("--out_prefix", required=True, help="Prefix for output files")
    args = parser.parse_args()

    adata = sc.read_h5ad(args.h5ad)
    coords = adata.obsm["spatial"]
    celltypes = adata.obs["dominant_celltype"]
    categories = celltypes.cat.categories

    # --- Plot 1: spatial scatter colored by dominant cell type ---
    fig, ax = plt.subplots(figsize=(7, 7))
    cmap = plt.get_cmap("tab10")
    for i, ct in enumerate(categories):
        mask = celltypes == ct
        ax.scatter(coords[mask, 0], coords[mask, 1], s=8, color=cmap(i % 10), label=ct, alpha=0.8)
    ax.invert_yaxis()
    ax.set_xlabel("Spatial X")
    ax.set_ylabel("Spatial Y")
    ax.set_title("Dominant cell type per spot")
    ax.legend(markerscale=2, fontsize=8, bbox_to_anchor=(1.02, 1), loc="upper left")
    ax.set_aspect("equal")
    fig.tight_layout()
    fig.savefig(f"{args.out_prefix}_spatial_celltypes.png", dpi=150)
    plt.close(fig)

    # --- Plot 2: neighborhood enrichment heatmap ---
    zscore = adata.uns["dominant_celltype_nhood_enrichment"]["zscore"]
    enrichment_df = pd.DataFrame(zscore, index=categories, columns=categories)

    fig, ax = plt.subplots(figsize=(7, 6))
    im = ax.imshow(enrichment_df.values, cmap="RdBu_r", vmin=-40, vmax=40)
    ax.set_xticks(range(len(categories)))
    ax.set_yticks(range(len(categories)))
    ax.set_xticklabels(categories, rotation=45, ha="right", fontsize=8)
    ax.set_yticklabels(categories, fontsize=8)
    ax.set_title("Neighborhood enrichment (z-score)")
    fig.colorbar(im, ax=ax, label="z-score")
    fig.tight_layout()
    fig.savefig(f"{args.out_prefix}_enrichment_heatmap.png", dpi=150)
    plt.close(fig)

    # --- Summary report ---
    with open(f"{args.out_prefix}_report.md", "w") as f:
        f.write("# Spatial transcriptomics pipeline results\n\n")
        f.write(f"Total spots analyzed: {adata.n_obs}\n\n")
        f.write("## Dominant cell type distribution\n\n")
        f.write(celltypes.value_counts().to_string())
        f.write("\n\n## Strongest co-localization pairs\n\n")
        off_diag = enrichment_df.values.copy()
        np.fill_diagonal(off_diag, -np.inf)
        top_idx = np.unravel_index(np.argmax(off_diag), off_diag.shape)
        f.write(f"{categories[top_idx[0]]} <-> {categories[top_idx[1]]}: z = {off_diag[top_idx]:.2f}\n")

    print(f"Visualization done: {args.out_prefix}_spatial_celltypes.png, {args.out_prefix}_enrichment_heatmap.png, {args.out_prefix}_report.md")

if __name__ == "__main__":
    main()
