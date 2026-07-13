process QC {
    conda "/home/arumuganainar.d/.conda/envs/spatial-tme"
    publishDir "${params.outdir}/qc", mode: 'copy'

    input:
    path h5ad_file

    output:
    path "filtered.h5ad", emit: filtered
    path "filtered_qc_summary.txt", emit: summary

    script:
    """
    run_qc.py --h5ad ${h5ad_file} --min_genes ${params.min_genes} --max_mito_pct ${params.max_mito_pct} --out filtered.h5ad
    """
}
