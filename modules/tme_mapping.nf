process TME_MAPPING {
    conda "/home/arumuganainar.d/.conda/envs/spatial-tme"
    publishDir "${params.outdir}/tme", mode: 'copy'

    input:
    path h5ad_file

    output:
    path "tme.h5ad", emit: tme_mapped
    path "tme_tme_summary.txt", emit: summary

    script:
    """
    run_tme_mapping.py --h5ad ${h5ad_file} --celltype_col celltype_major --n_perms 1000 --out tme.h5ad
    """
}
