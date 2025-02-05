process QIIME2_TAXA_BAR_PLOT {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::qiime2=2024.2"
    container "quay.io/qiime2/core:2024.2"

    input:
    tuple val(meta), path(taxonomy)
    path table
    path metadata

    output:
    tuple val(meta), path("*-taxa-bar-plots.qzv"), emit: barplot
    path "versions.yml", emit: versions

    script:
    """
    qiime taxa barplot \\
        --i-table $table \\
        --i-taxonomy $taxonomy \\
        --m-metadata-file $metadata \\
        --o-visualization ${meta.id}-taxa-bar-plots.qzv

    echo 'QIIME2_VERSION: ' \$(qiime --version) > versions.yml
    """
}
