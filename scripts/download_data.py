#!/usr/bin/env python
"""
Fetch the two pipeline inputs into data/ (run once, from the spatial conda env,
on a node with internet access):

  data/visium_breast_cancer.h5ad  10x Visium V1_Breast_Cancer_Block_A_Section_1 (via Squidpy)
  data/wu2021_breast_atlas.h5ad   Wu et al. 2021 breast cancer atlas, 100,064 cells (CELLxGENE)

The atlas is resolved through the CELLxGENE curation API by its stable dataset ID,
so a re-published version is picked up automatically; the exact file URL is printed.
"""
import argparse
import json
import shutil
import urllib.request
from pathlib import Path

VISIUM_SAMPLE = "V1_Breast_Cancer_Block_A_Section_1"
CELLXGENE_COLLECTION = "dea97145-f712-431c-a223-6b5f565f362a"  # doi:10.1038/s41588-021-00911-1
CELLXGENE_DATASET = "9fddb063-056d-4202-8b8a-4b0ee531d3ce"     # full atlas, all major cell types
API = "https://api.cellxgene.cziscience.com/curation/v1"


def download_visium(out):
    import squidpy as sq

    adata = sq.datasets.visium(VISIUM_SAMPLE, base_dir=out.parent / "visium")
    adata.var_names_make_unique()
    adata.write_h5ad(out)
    print(f"Visium: {adata.n_obs} spots x {adata.n_vars} genes -> {out}")


def download_reference(out):
    url = f"{API}/collections/{CELLXGENE_COLLECTION}/datasets/{CELLXGENE_DATASET}"
    with urllib.request.urlopen(url) as r:
        dataset = json.load(r)
    h5ad_url = next(a["url"] for a in dataset["assets"] if a["filetype"] == "H5AD")
    print(f"Reference: {dataset['title']} ({dataset['cell_count']} cells)\n  {h5ad_url}")

    tmp = out.with_suffix(".part")
    with urllib.request.urlopen(h5ad_url) as r, open(tmp, "wb") as f:
        shutil.copyfileobj(r, f, length=1 << 20)
    tmp.rename(out)
    print(f"  -> {out}")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--outdir", default=Path(__file__).resolve().parent.parent / "data", type=Path)
    parser.add_argument("--skip_reference", action="store_true", help="only fetch the Visium sample (~850 MB saved)")
    args = parser.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    visium = args.outdir / "visium_breast_cancer.h5ad"
    if visium.exists():
        print(f"Visium: {visium} exists, skipping")
    else:
        download_visium(visium)

    reference = args.outdir / "wu2021_breast_atlas.h5ad"
    if args.skip_reference:
        pass
    elif reference.exists():
        print(f"Reference: {reference} exists, skipping")
    else:
        download_reference(reference)


if __name__ == "__main__":
    main()
