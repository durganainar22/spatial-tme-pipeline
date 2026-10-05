process SVG_DETECTION {
    label 'process_medium'
    conda "${moduleDir}/../envs/spatial.yml"
    publishDir "${params.outdir}/svg", mode: 'copy'

    input:
    path h5ad_file

    output:
    path "svg.h5ad", emit: svg_flagged
    path "svg_svg_summary.txt", emit: summary

    script:
    """
    run_svg_detection.py --h5ad ${h5ad_file} --n_top_svgs ${params.n_top_svgs} --out svg.h5ad
    """

    stub:
    """
    touch svg.h5ad svg_svg_summary.txt
    """
}
