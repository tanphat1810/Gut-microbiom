process BETA_EMPEROR {
    tag "qiime2_beta_emperor"

    conda "bioconda::qiime2=2023.7"
    container "quay.io/qiime2/amplicon:2024.5"

    input:
    path pcoa_results  // Kết quả từ BETA_PCOA
    path metadata_file // Metadata để gán nhóm mẫu

    output:
    path "*_emperor.qzv", emit: emperor_plots
    publishDir "qiime_out", mode: 'copy'

    script:
    """
    # 🔹 Tạo biểu đồ PCoA bằng Emperor
    for pcoa in ${pcoa_results}; do
        basename=\$(basename \$pcoa .qza)
        qiime emperor plot \\
            --i-pcoa \$pcoa \\
            --m-metadata-file ${metadata_file} \\
            --o-visualization \${basename}_emperor.qzv
    done
    """
}

