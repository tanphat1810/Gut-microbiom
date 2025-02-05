process QIIME2_CHIMERA_FREE_REP_SEQS_SUMMARY {
    tag "$meta.id"
    label 'process_low'

    conda "bioconda::qiime2=2024.2"
    container "quay.io/qiime2/core:2024.2"

    input:
    tuple val(meta), path(chimera_free_rep_seqs)

    output:
    tuple val(meta), path("*-rep-seqs.qzv"), emit: rep_seqs_summary_qzv
    path "versions.yml", emit: versions

    script:
    """
    qiime feature-table tabulate-seqs \\
        --i-data $chimera_free_rep_seqs \\
        --o-visualization ${meta.id}-rep-seqs.qzv

    echo 'QIIME2_VERSION: ' \$(qiime --version) > versions.yml
    """
}
