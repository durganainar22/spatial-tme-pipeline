process REFERENCE_SIGNATURES {
    label 'gpu'
    conda "${moduleDir}/../envs/cell2location.yml"
    publishDir "${params.outdir}/reference", mode: 'copy'

    input:
    path ref_h5ad

    output:
    path "reference_signatures.h5ad", emit: signatures
    path "reference_model", emit: model

    script:
    """
    estimate_reference_signatures.py --ref_h5ad ${ref_h5ad} --celltype_col ${params.celltype_col} --batch_col ${params.batch_col} --max_epochs ${params.ref_max_epochs} --accelerator ${params.accelerator} --out .
    """

    stub:
    """
    mkdir reference_model
    touch reference_signatures.h5ad
    """
}
