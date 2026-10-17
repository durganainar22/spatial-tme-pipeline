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

Each stage is an independent Nextflow process, chained via explicit channel outputs (Process.out.label). The pipeline runs across two isolated conda environments -- a lightweight scanpy/squidpy environment for stages 1-3 and 5-6, and a GPU-enabled cell2location environment for stage 4 -- defined in `envs/` and selected per process via Nextflow's conda directive.

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
    |-- nextflow.config            # params, profiles (conda, mamba, slurm)
    |-- conf/
    |   |-- base.config            # per-label cpus/memory/time, capped by --max_*
    |   `-- slurm.config           # SLURM executor + GPU request for cell2location steps
    |-- envs/
    |   |-- spatial.yml            # scanpy/squidpy env (pinned)
    |   `-- cell2location.yml      # cell2location/scvi-tools/torch env (pinned)
    |-- modules/                   # one .nf file per pipeline stage
    |   |-- qc.nf
    |   |-- normalize.nf
    |   |-- spatially_variable_genes.nf
    |   |-- reference_signatures.nf
    |   |-- deconvolution.nf
    |   |-- tme_mapping.nf
    |   `-- visualize.nf
    |-- bin/                       # Python scripts called by each process
    |-- scripts/download_data.py   # fetches the Visium sample + reference atlas into data/
    |-- data/                      # inputs (gitignored)
    `-- results/                   # pipeline outputs (gitignored)

## Running the pipeline

Requires Nextflow >= 24.04 and conda/mamba. Each process builds its environment from `envs/` on first use.

1. Fetch the inputs (once, on a machine with internet access; the atlas is ~850 MB):

        conda env create -f envs/spatial.yml
        conda activate spatial-tme
        python scripts/download_data.py

2. Run on a SLURM cluster (GPU needed for the two cell2location steps):

        nextflow run main.nf -profile conda,slurm

   Queues default to `short` / `gpu`; change them with `--cpu_queue`, `--gpu_queue`, `--gpu_type v100`.

   Or run locally with a GPU: `nextflow run main.nf -profile conda --max_memory 16.GB`

The first run trains the reference signatures (REFERENCE_SIGNATURES) and publishes them to `results/reference/`. Later runs can skip that step:

    nextflow run main.nf -profile conda,slurm --ref_signatures results/reference/reference_signatures.h5ad

Add `-resume` to reuse cached steps. `-stub` runs the whole workflow with placeholder outputs to check the wiring without data or a GPU.

## Author

Durga Gomathi Arumuganainar -- MS Bioinformatics, Northeastern University
