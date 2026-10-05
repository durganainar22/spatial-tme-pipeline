process DECONVOLUTION {
    label 'gpu'
    conda "${moduleDir}/../envs/cell2location.yml"
    publishDir "${params.outdir}/deconvolution", mode: 'copy'

    input:
    path spatial_h5ad
    path ref_signatures

    output:
    path "deconvolution.h5ad", emit: deconvolved

    script:
    """
    run_deconvolution.py --spatial_h5ad ${spatial_h5ad} --ref_signatures ${ref_signatures} --celltype_col ${params.celltype_col} --max_epochs ${params.max_epochs} --accelerator ${params.accelerator} --out deconvolution.h5ad
    """

    stub:
    """
    touch deconvolution.h5ad
    """
}
