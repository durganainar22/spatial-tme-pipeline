process TME_MAPPING {
    label 'process_medium'
    conda "${moduleDir}/../envs/spatial.yml"
    publishDir "${params.outdir}/tme", mode: 'copy'

    input:
    path h5ad_file

    output:
    path "tme.h5ad", emit: tme_mapped
    path "tme_tme_summary.txt", emit: summary

    script:
    """
    run_tme_mapping.py --h5ad ${h5ad_file} --n_perms ${params.n_perms} --out tme.h5ad
    """

    stub:
    """
    touch tme.h5ad tme_tme_summary.txt
    """
}
