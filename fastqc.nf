process FASTQC {
    tag "$meta.id"
    label 'process_low'

    conda "bioconda::fastqc=0.12.1"
    container "quay.io/biocontainers/fastqc:0.12.1--hdfd78af_0"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.html"), emit: html
    tuple val(meta), path("*.zip"), emit: zip
    path "versions.yml", emit: versions

    script:
    """
    fastqc --threads $task.cpus --outdir . $reads
    echo 'FASTQC_VERSION: ' \$(fastqc --version) > versions.yml
    """
}
