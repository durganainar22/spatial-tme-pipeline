process NORMALIZE {
    conda "/home/arumuganainar.d/.conda/envs/spatial-tme"
    publishDir "${params.outdir}/normalize", mode: 'copy'

    input:
    path h5ad_file

    output:
    path "normalized.h5ad", emit: normalized
    path "normalized_norm_summary.txt", emit: summary

    script:
    """
    run_normalization.py --h5ad ${h5ad_file} --n_top_genes ${params.n_top_genes} --out normalized.h5ad
    """
}
