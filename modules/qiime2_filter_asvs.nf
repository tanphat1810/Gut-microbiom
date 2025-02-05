process QIIME2_FILTER_ASVS {
    tag "$meta.id"
    label 'process_low'

    conda "bioconda::qiime2=2023.9"
    container "quay.io/qiime2/core:2023.9"

    input:
    tuple val(meta), path(table_qza)

    output:
    tuple val(meta), path("filtered-table.qza"), emit: filtered_table
    path "versions.yml", emit: versions

    script:
    """
    qiime feature-table filter-features \\
      --i-table $table_qza \\
      --p-min-frequency ${params.min_frequency} \\
      --p-min-samples ${params.min_samples} \\
      --o-filtered-table filtered-table.qza

    echo 'QIIME2_VERSION: ' \$(qiime --version) > versions.yml
    """
}
