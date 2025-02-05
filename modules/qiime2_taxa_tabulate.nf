process QIIME2_TAXA_TABULATE {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::qiime2=2024.2"
    container "quay.io/qiime2/core:2024.2"

    input:
    tuple val(meta), path(taxonomy)

    output:
    tuple val(meta), path("*-taxonomy-tabulate.qzv"), emit: tabulate
    path "versions.yml", emit: versions

    script:
    """
    qiime metadata tabulate \\
        --m-input-file $taxonomy \\
        --o-visualization ${meta.id}-taxonomy-tabulate.qzv

    echo 'QIIME2_VERSION: ' \$(qiime --version) > versions.yml
    """
}
