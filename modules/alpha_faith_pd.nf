process FAITH_PD {
    tag "qiime2_faith_pd"

    conda "bioconda::qiime2=2023.7"
    container "quay.io/qiime2/amplicon:2024.5"

    input:
    path filtered_sample_table   // Bảng feature table
    path rooted_tree     // Cây phát sinh loài

    output:
    path "faith_pd.qza", emit: faith_pd
    publishDir "qiime_out", mode: 'copy'

    script:
    """
    qiime diversity alpha-phylogenetic \\
        --i-table ${filtered_sample_table} \\
        --i-phylogeny ${rooted_tree} \\
        --p-metric faith_pd \\
        --o-alpha-diversity faith_pd.qza
    """
}

