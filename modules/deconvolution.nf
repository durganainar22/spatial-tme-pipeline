process DECONVOLUTION {
    conda "/home/arumuganainar.d/.conda/envs/cell2location"
    publishDir "${params.outdir}/deconvolution", mode: 'copy'

    input:
    path spatial_h5ad
    path ref_signatures

    output:
    path "deconvolution.h5ad", emit: deconvolved

    script:
    """
    run_deconvolution.py --spatial_h5ad ${spatial_h5ad} --ref_signatures ${ref_signatures} --celltype_col celltype_major --max_epochs 500 --out deconvolution.h5ad
    """
}
