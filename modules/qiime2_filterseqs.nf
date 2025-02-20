process FILTER_SEQS {
    tag "qiime2_filter_seqs"

    conda "bioconda::qiime2=2023.7"
    container "quay.io/qiime2/amplicon:2024.5"

    input:
    path rep_seqs          
    path filtered_table      

    output:
    path "filtered_rep_seqs.qza", emit: filtered_representative_seqs

    script:
    """
    qiime feature-table filter-seqs \\
        --i-data ${rep_seqs} \\
        --i-table ${filtered_table} \\
        --o-filtered-data filtered_rep_seqs.qza
    """
}

