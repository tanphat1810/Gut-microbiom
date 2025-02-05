process ALPHA_DIVERSITY_VISUALIZATION {
    tag "$meta.id"
    label 'process_low'

    conda "bioconda::qiime2=2023.5"
    container "quay.io/biocontainers/qiime2:2023.5"

    input:
    tuple val(meta), path(alpha_diversity), path(metadata)

    output:
    tuple val(meta), path("*_alpha_significance.qzv"), emit: alpha_significance
    path "versions.yml", emit: versions

    script:
    """
    #  Tạo báo cáo thống kê sự khác biệt alpha diversity giữa nhóm
    qiime diversity alpha-group-significance \
        --i-alpha-diversity $alpha_diversity \
        --m-metadata-file $metadata \
        --o-visualization alpha_significance.qzv

    #  Lưu phiên bản phần mềm
    echo 'QIIME2_VERSION: ' \$(qiime --version) > versions.yml
    """
}
