process BETA_DIVERSITY_VISUALIZATION {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::qiime2=2023.5"
    container "quay.io/biocontainers/qiime2:2023.5"

    input:
    tuple val(meta), path(beta_distance_matrix), path(metadata)

    output:
    tuple val(meta), path("*_emperor.qzv"), emit: beta_emperor_plot
    path "versions.yml", emit: versions

    script:
    """
    #  Tạo biểu đồ PCoA bằng Emperor
    qiime emperor plot \
        --i-pcoa $beta_distance_matrix \
        --m-metadata-file $metadata \
        --o-visualization beta_emperor.qzv

    #  Lưu phiên bản phần mềm
    echo 'QIIME2_VERSION: ' \$(qiime --version) > versions.yml
    """
}
