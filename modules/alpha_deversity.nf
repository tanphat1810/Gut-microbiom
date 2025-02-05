process ALPHA_DIVERSITY {
    tag "$meta.id"
    label 'process_low'

    conda "bioconda::qiime2=2023.5"
    container "quay.io/biocontainers/qiime2:2023.5"

    input:
    tuple val(meta), path(table)

    output:
    tuple val(meta), path("*_alpha.qza"), emit: alpha_diversity
    path "versions.yml", emit: versions

    script:
    """
    # Tính toán các chỉ số đa dạng alpha phổ biến
    qiime diversity alpha \
        --i-table $table \
        --p-metric shannon \
        --o-alpha-diversity shannon_alpha.qza

    qiime diversity alpha \
        --i-table $table \
        --p-metric simpson \
        --o-alpha-diversity simpson_alpha.qza

    qiime diversity alpha \
        --i-table $table \
        --p-metric chao1 \
        --o-alpha-diversity chao1_alpha.qza

    qiime diversity alpha \
        --i-table $table \
        --p-metric pielou_e \
        --o-alpha-diversity pielou_evenness_alpha.qza

    qiime diversity alpha \
        --i-table $table \
        --p-metric faith_pd \
        --o-alpha-diversity faith_pd_alpha.qza

    # Lưu phiên bản phần mềm
    echo 'QIIME2_VERSION: ' \$(qiime --version) > versions.yml
    """
}
