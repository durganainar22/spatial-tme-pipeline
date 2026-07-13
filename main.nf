#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.h5ad = "${projectDir}/test_data/visium_sample.h5ad"
params.outdir = "${projectDir}/results"
params.min_genes = 200
params.max_mito_pct = 20.0
params.n_top_genes = 2000
params.n_top_svgs = 50
params.ref_signatures = "${projectDir}/results/reference_model/reference_signatures.h5ad"

include { QC } from './modules/qc.nf'
include { NORMALIZE } from './modules/normalize.nf'
include { SVG_DETECTION } from './modules/spatially_variable_genes.nf'
include { DECONVOLUTION } from './modules/deconvolution.nf'
include { TME_MAPPING } from './modules/tme_mapping.nf'
include { VISUALIZE } from './modules/visualize.nf'

workflow {
    h5ad_ch = Channel.fromPath(params.h5ad)
    ref_ch = Channel.fromPath(params.ref_signatures)

    QC(h5ad_ch)
    NORMALIZE(QC.out.filtered)
    SVG_DETECTION(NORMALIZE.out.normalized)
    DECONVOLUTION(QC.out.filtered, ref_ch)
    TME_MAPPING(DECONVOLUTION.out.deconvolved)
    VISUALIZE(TME_MAPPING.out.tme_mapped)
}
