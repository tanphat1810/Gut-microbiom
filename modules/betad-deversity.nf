process BETA_DIVERSITY {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::qiime2=2023.5"
    container "quay.io/biocontainers/qiime2:2023.5"

    input:
    tuple val(meta), path(table), path(phylogeny)

    output:
    tuple val(meta), path("*_beta.qza"), emit: beta_diversity
    path "versions.yml", emit: versions

    script:
    """
    # ✅ Tính toán chỉ số beta diversity không dựa trên phát sinh loài
    qiime diversity beta \
        --i-table $table \
        --p-metric braycurtis \
        --o-distance-matrix braycurtis_beta.qza

    qiime diversity beta \
        --i-table $table \
        --p-metric jaccard \
        --o-distance-matrix jaccard_beta.qza

    # ✅ Tính toán chỉ số beta diversity dựa trên cây phát sinh loài
    qiime diversity beta-phylogenetic \
        --i-table $table \
        --i-phylogeny $phylogeny \
        --p-metric unweighted_unifrac \
        --o-distance-matrix unweighted_unifrac_beta.qza

    qiime diversity beta-phylogenetic \
        --i-table $table \
        --i-phylogeny $phylogeny \
        --p-metric weighted_unifrac \
        --o-distance-matrix weighted_unifrac_beta.qza

    # ✅ Lưu phiên bản phần mềm
    echo 'QIIME2_VERSION: ' \$(qiime --version) > versions.yml
    """
}
