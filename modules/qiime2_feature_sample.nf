process FILTER_SAMPLES {
    tag "qiime2_filter_samples"

    conda "bioconda::qiime2=2023.7"
    container "quay.io/qiime2/amplicon:2024.5"

    input:
    path filtered_feature_table  // Nhận đầu vào là bảng feature table đã lọc

    output:
    path "filtered_sample_table.qza", emit: filtered_sample_table
    publishDir "qiime_out", mode: 'copy'
    script:
    """
    qiime feature-table filter-samples \\
        --i-table ${filtered_feature_table} \\
        --p-min-frequency ${params.filtersample.min_sample_freq} \\
        --o-filtered-table filtered_sample_table.qza
    """
}

