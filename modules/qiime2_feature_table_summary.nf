process QIIME2_FEATURE_TABLE_SUMMARY {
    tag "$meta.id"
    label 'process_low'

    conda "bioconda::qiime2=2023.9"
    container "quay.io/qiime2/core:2023.9"

    input:
    tuple val(meta), path(table_qza)

    output:
    tuple val(meta), path("table-summary.qzv"), emit: table_summary_qzv
    path "versions.yml", emit: versions

    script:
    """
    qiime feature-table summarize \\
      --i-table $table_qza \\
      --o-visualization table-summary.qzv

    echo 'QIIME2_VERSION: ' \$(qiime --version) > versions.yml
    """
}
