process QIIME2_FILTER_CHIMERIC_ASVS {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::qiime2=2024.2"
    container "quay.io/qiime2/core:2024.2"

    input:
    tuple val(meta), path(feature_table), path(rep_seqs)

    output:
    tuple val(meta), path("*-filtered-table.qza"), emit: filtered_table
    tuple val(meta), path("*-filtered-rep-seqs.qza"), emit: filtered_rep_seqs
    path "versions.yml", emit: versions

    script:
    """
    qiime chimera.vsearch \\
        --i-table $feature_table \\
        --i-sequences $rep_seqs \\
        --o-chimeric-table ${meta.id}-chimeric-table.qza \\
        --o-nonchimeric-table ${meta.id}-filtered-table.qza \\
        --o-chimeric-sequences ${meta.id}-chimeric-rep-seqs.qza \\
        --o-nonchimeric-sequences ${meta.id}-filtered-rep-seqs.qza

    echo 'QIIME2_VERSION: ' \$(qiime --version) > versions.yml
    """
}
