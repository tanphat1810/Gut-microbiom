process ALPHA_METRICS {
    tag "qiime2_alpha_metrics"

    conda "bioconda::qiime2=2023.7"
    container "quay.io/qiime2/amplicon:2024.5"

    input:
    path feature_table   // Bảng feature table

    output:
    path "observed_otus.qza", emit: observed_otus
    path "chao1.qza", emit: chao1
    path "shannon.qza", emit: shannon
    path "simpson.qza", emit: simpson
    path "pielou_e.qza", emit: pielou
    publishDir "qiime_out", mode: 'copy'

    script:
    """
    for metric in observed_otus chao1 shannon simpson pielou_e; do
        qiime diversity alpha \\
            --i-table ${feature_table} \\
            --p-metric \${metric} \\
            --o-alpha-diversity \${metric}.qza
    done
    """
}

