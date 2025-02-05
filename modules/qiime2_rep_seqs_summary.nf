process QIIME2_REP_SEQS_SUMMARY {
    tag "$meta.id"
    label 'process_low'

    conda "bioconda::qiime2=2023.9"
    container "quay.io/qiime2/core:2023.9"

    input:
    tuple val(meta), path(rep_seqs_qza)

    output:
    tuple val(meta), path("rep-seqs-summary.qzv"), emit: rep_seqs_summary_qzv
    path "versions.yml", emit: versions

    script:
    """
    qiime feature-table tabulate-seqs \\
      --i-data $rep_seqs_qza \\
      --o-visualization rep-seqs-summary.qzv

    echo 'QIIME2_VERSION: ' \$(qiime --version) > versions.yml
    """
}
