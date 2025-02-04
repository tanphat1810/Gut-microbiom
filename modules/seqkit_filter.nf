process SEQKIT_FILTER {
    tag "$meta.id"
    label 'process_low'

    conda "bioconda::seqkit=2.4.0"
    container "quay.io/biocontainers/seqkit:2.4.0--h9ee0642_0"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("filtered/*.fastq.gz"), emit: filtered_reads
    path "versions.yml", emit: versions

    script:
    """
    mkdir -p filtered
    seqkit fq -q 20 -o filtered/${meta.id}.filtered.fastq.gz $reads

    echo 'SEQKIT_VERSION: ' \$(seqkit version | grep -oE '[0-9.]+') > versions.yml
    """
}
