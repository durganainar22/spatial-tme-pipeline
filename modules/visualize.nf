process VISUALIZE {
    conda "/home/arumuganainar.d/.conda/envs/spatial-tme"
    publishDir "${params.outdir}/visualize", mode: 'copy'

    input:
    path h5ad_file

    output:
    path "pipeline_spatial_celltypes.png"
    path "pipeline_enrichment_heatmap.png"
    path "pipeline_report.md"

    script:
    """
    run_visualization.py --h5ad ${h5ad_file} --out_prefix pipeline
    """
}
