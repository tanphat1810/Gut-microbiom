process FILTER_FEATURES {
    tag "filter_feature_sample"

    conda "bioconda::qiime2=2023.7"
    container "quay.io/qiime2/amplicon:2024.5"
    input:
    path feature_table    
    output:
    path "filtered_features_table.qza", emit: filtered_feature_table
    publishDir "qiime_out", mode: 'copy'
    script:
    """
    qiime feature-table filter-features \\
        --i-table ${feature_table} \\
        --p-min-frequency ${params.filter.min_feature_freq} \\
        --o-filtered-table filtered_features_table.qza  
    """
}

