process QIIME2_CHIMERA_FREE_TABLE_SUMMARY {
    tag "$meta.id"
    label 'process_low'

    conda "bioconda::qiime2=2024.2"
    container "quay.io/qiime2/core:2024.2"

    input:
    tuple val(meta), path(chimera_free_table)

    output:
    tuple val(meta), path("*-summary.qzv"), emit: table_summary_qzv
    path "versions.yml", emit: versions

    script:
    """
    qiime feature-table summarize \\
        --i-table $chimera_free_table \\
        --o-visualization ${meta.id}-summary.qzv

    echo 'QIIME2_VERSION: ' \$(qiime --version) > versions.yml
    """
}
