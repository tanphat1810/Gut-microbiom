process SEQKIT_FILTER {
    tag "$meta.id"
    label 'process_low'

    conda "bioconda::seqkit=2.3.1"
    container "quay.io/biocontainers/seqkit:2.3.1--h9ee0642_0"

    input:
    tuple val(meta), path(read)

    output:
    tuple val(meta), path("${meta.id}.filtered.fastq.gz"), emit: filtered
    tuple val(meta), path("${meta.id}.log"), emit: log
    publishDir "seqkit_result", mode: 'copy'
   
    script:
    """
    seqkit seq -Q 20 -j 4 "${read}" | gzip > "${meta.id}.filtered.fastq.gz"

    echo "SeqKit filtering completed for ${meta.id}" > "${meta.id}.log"
    """
} 
