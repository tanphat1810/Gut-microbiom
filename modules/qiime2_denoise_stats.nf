process QIIME2_DENOISE_STATS {
    tag "$meta.id"
    label 'process_low'

    conda "bioconda::qiime2=2023.9"
    container "quay.io/qiime2/core:2023.9"

    input:
    tuple val(meta), path(stats_qza)

    output:
    tuple val(meta), path("stats.qzv"), emit: stats_qzv
    path "versions.yml", emit: versions

    script:
    """
    qiime metadata tabulate \\
      --m-input-file $stats_qza \\
      --o-visualization stats.qzv

    echo 'QIIME2_VERSION: ' \$(qiime --version) > versions.yml
    """
}
