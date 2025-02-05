process QIIME2_TAXONOMIC_CLASSIFICATION {
    tag "$meta.id"
    label 'process_high'

    conda "bioconda::qiime2=2024.2"
    container "quay.io/qiime2/core:2024.2"

    input:
    tuple val(meta), path(rep_seqs)
    path classifier
    path reference_database

    output:
    tuple val(meta), path("*-taxonomy.qza"), emit: taxonomy
    tuple val(meta), path("*-taxonomy.tsv"), emit: taxonomy_tsv
    path "versions.yml", emit: versions

    script:
    """
    qiime feature-classifier classify-sklearn \\
        --i-classifier $classifier \\
        --i-reads $rep_seqs \\
        --o-classification ${meta.id}-taxonomy.qza

    qiime tools export \\
        --input-path ${meta.id}-taxonomy.qza \\
        --output-path taxonomy_export

    mv taxonomy_export/taxonomy.tsv ${meta.id}-taxonomy.tsv

    echo 'QIIME2_VERSION: ' \$(qiime --version) > versions.yml
    """
}
