# Spatial transcriptomics pipeline for tumor microenvironment analysis

An end-to-end Nextflow DSL2 pipeline that processes 10x Visium spatial transcriptomics data, deconvolves cell type composition per spot using a public single-cell reference atlas, and characterizes spatial co-localization patterns within the tumor microenvironment.

## Pipeline overview

    Raw Visium data (.h5ad)
            |
            v
       QC ---------------> filter spots by gene count, mitochondrial %
            |
            v
       NORMALIZE ---------> log-normalize, select highly variable genes
            |
            v
       SVG_DETECTION ------> Moran's I spatial autocorrelation (Squidpy)
            |
            v
       DECONVOLUTION ------> per-spot cell type proportions (cell2location, GPU)
            |                reference: Wu et al. 2021, Nature Genetics (breast cancer atlas)
            v
       TME_MAPPING --------> neighborhood enrichment between cell types
            |
            v
       VISUALIZE ----------> spatial plots, enrichment heatmap, summary report

Each stage is an independent Nextflow process, chained via explicit channel outputs (Process.out.label). The pipeline runs across two isolated conda environments -- a lightweight scanpy/squidpy environment for stages 1-3 and 5-6, and a GPU-enabled cell2location environment for stage 4 -- selected automatically per-process via Nextflow's conda directive.

## Key results

On the demo dataset (10x Visium breast cancer, 3798 spots):

- **Spatially variable genes** cleanly separated into three biological compartments: stromal (COL1A1, COL1A2, DCN), immune (immunoglobulin genes), and epithelial (keratins) -- all with FDR-corrected p ~ 0.
- **Deconvolution** estimated proportions across 9 cell types (Cancer Epithelial, Normal Epithelial, CAFs, T-cells, Myeloid, Endothelial, PVL, B-cells, Plasmablasts) using a reference-trained regression model.
- **Neighborhood enrichment** confirmed known tissue architecture without being told what to look for: Cancer Epithelial forms a segregated compartment (strongly negative co-localization with every other cell type), while Endothelial, PVL, and CAFs cluster together (vascularized stroma), and Myeloid/T-cells/B-cells cluster together (immune infiltrate).

## Tech stack

| Component | Tool |
|---|---|
| Workflow orchestration | Nextflow (DSL2) |
| Spatial analysis | Squidpy |
| Single-cell analysis | Scanpy |
| Cell type deconvolution | Cell2location (PyTorch, scvi-tools) |
| Compute | Northeastern Discovery HPC (SLURM), NVIDIA V100 GPU |

## Data

- **Spatial data**: 10x Genomics Visium, breast cancer FFPE section (V1_Breast_Cancer_Block_A_Section_1, via Squidpy's built-in dataset loader). Used as a placeholder demo dataset.
- **Reference atlas**: Wu et al. 2021, Nature Genetics -- single-cell and spatially resolved atlas of human breast cancers, 100,064 cells across 9 major cell types, accessed via CZ CELLxGENE Discover (https://cellxgene.cziscience.com/).

## Repository structure

    .
    |-- main.nf                    # workflow definition, wires all processes together
    |-- nextflow.config            # conda.enabled + global config
    |-- modules/                   # one .nf file per pipeline stage
    |   |-- qc.nf
    |   |-- normalize.nf
    |   |-- spatially_variable_genes.nf
    |   |-- deconvolution.nf
    |   |-- tme_mapping.nf
    |   `-- visualize.nf
    |-- bin/                       # Python scripts called by each process
    |   |-- run_qc.py
    |   |-- run_normalization.py
    |   |-- run_svg_detection.py
    |   |-- estimate_reference_signatures.py
    |   |-- run_deconvolution.py
    |   |-- run_tme_mapping.py
    |   `-- run_visualization.py
    |-- test_data/                 # input spatial data + reference atlas
    `-- results/                   # pipeline outputs (gitignored)

## Running the pipeline

Requires two conda environments: spatial-tme (scanpy, squidpy, nextflow) and cell2location (torch, cell2location, scvi-tools). See bin/ scripts for exact package requirements.

    nextflow run main.nf

Reference signatures must be pre-trained once via bin/estimate_reference_signatures.py before the first full pipeline run (see results/reference_model/ for the expected output location).

## Author

Durga Gomathi Arumuganainar -- MS Bioinformatics, Northeastern University
