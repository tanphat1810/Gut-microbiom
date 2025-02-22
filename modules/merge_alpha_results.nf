process MERGE_ALPHA_RESULTS {
    tag "qiime2_merge_alpha_diversity"

    conda "bioconda::qiime2=2023.7"
    container "quay.io/qiime2/amplicon:2024.5"

    input:
    path observed_otus
    path chao1
    path shannon
    path simpson
    path pielou
    path faith_pd

    output:
    path "alpha_diversity_metrics.qzv", emit: alpha_metrics
    publishDir "qiime_out", mode: 'copy'

    script:
    """
    qiime metadata tabulate \\
        --m-input-file ${observed_otus} \\
        --m-input-file ${chao1} \\
        --m-input-file ${shannon} \\
        --m-input-file ${simpson} \\
        --m-input-file ${pielou} \\
        --m-input-file ${faith_pd} \\
        --o-visualization alpha_diversity_metrics.qzv
    """
}

