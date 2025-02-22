process CLASSIFY_TAXONOMY {
    tag "qiime2_classify_sklearn"

    conda "bioconda::qiime2=2023.7"
    container "quay.io/qiime2/amplicon:2024.5"

    input:
    path filtered_rep_seqs  // Trình tự đại diện đã được lọc
    path classifier         // Bộ phân loại Silva

    output:
    path "taxonomy.qza", emit: classified_taxonomy
    path "taxonomy.qzv", emit: classified_taxonomy_viz
    publishDir "qiime_out", mode: 'copy'

    script:
    """
    qiime feature-classifier classify-sklearn \\
        --i-classifier $classifier \\
        --i-reads $filtered_rep_seqs \\
        --o-classification taxonomy.qza

    qiime metadata tabulate \\
        --m-input-file taxonomy.qza \\
        --o-visualization taxonomy.qzv
    """
}

