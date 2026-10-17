#!/usr/bin/env nextflow

include { QC } from './modules/qc.nf'
include { NORMALIZE } from './modules/normalize.nf'
include { SVG_DETECTION } from './modules/spatially_variable_genes.nf'
include { REFERENCE_SIGNATURES } from './modules/reference_signatures.nf'
include { DECONVOLUTION } from './modules/deconvolution.nf'
include { TME_MAPPING } from './modules/tme_mapping.nf'
include { VISUALIZE } from './modules/visualize.nf'

workflow {
    h5ad_ch = Channel.fromPath(params.h5ad, checkIfExists: true)

    // reuse trained reference signatures if given, otherwise train them from the atlas
    if (params.ref_signatures) {
        ref_ch = Channel.fromPath(params.ref_signatures, checkIfExists: true).first()
    } else {
        REFERENCE_SIGNATURES(Channel.fromPath(params.ref_h5ad, checkIfExists: true))
        ref_ch = REFERENCE_SIGNATURES.out.signatures.first()
    }

    QC(h5ad_ch)
    NORMALIZE(QC.out.filtered)
    SVG_DETECTION(NORMALIZE.out.normalized)
    DECONVOLUTION(QC.out.filtered, ref_ch)
    TME_MAPPING(DECONVOLUTION.out.deconvolved)
    VISUALIZE(TME_MAPPING.out.tme_mapped)
}
