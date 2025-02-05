process QIIME2_IMPORT {
    tag "$meta.id"
    label 'process_low'

    conda "bioconda::qiime2=2023.9"
    container "quay.io/biocontainers/qiime2:2023.9--py38_0"

    input:
    tuple val(meta), path(metadata)

    output:
    tuple val(meta), path("demux.qza"), emit: qiime2_artifact

    script:
    """
    qiime tools import \\
        --type 'SampleData[Sequences]' \\
        --input-path $metadata \\
        --output-path demux.qza
    """
}
