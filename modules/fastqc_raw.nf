process FASTQC_RAW {
    tag "$meta.id"
    label 'process_low'

    conda "bioconda::fastqc=0.12.1"
    container "quay.io/biocontainers/fastqc:0.12.1--hdfd78af_0"

    input:
    tuple val(meta), path(sample)

    output:
    tuple val(meta), path("*fastqc.html"), emit: html
    tuple val(meta), path("*fastqc.zip"), emit: zip
    publishDir "qiime_out", mode: 'copy'
    script:
    """
    fastqc --threads 4 $sample
    """
}

