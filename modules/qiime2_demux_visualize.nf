process QIIME2_DEMUX_VISUALIZE {
    tag "$meta.id"
    label 'process_low'

    conda "bioconda::qiime2=2023.7"
    container "quay.io/biocontainers/qiime2-core:2023.7"

    input:
    tuple val(meta), path(demux_qza)

    output:
    tuple val(meta), path("*.qzv"), emit: demux_qzv
    path "versions.yml", emit: versions

    script:
    """
    qiime demux summarize \\
        --i-data $demux_qza \\
        --o-visualization demux.qzv

    echo 'QIIME2_VERSION: ' \$(qiime --version) > versions.yml
    """
}
