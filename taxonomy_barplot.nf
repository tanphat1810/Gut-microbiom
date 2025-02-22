process BARPLOT {
    tag "qiime2_barplot"

    conda "bioconda::qiime2=2024.5"
    container 'quay.io/qiime2/amplicon:2024.5'  

    input:
    path taxonomy_qza       // Nhận đầu vào từ CLASSIFY_TAXONOMY
    path feature_sample     // Nhận đầu vào từ FILTER_SAMPLES
    path metadata           // Nhận metadata từ nextflow.config

    output:
    path "taxonomy_barplot.qzv", emit: barplot_qzv  
    publishDir "qiime_out", mode: 'copy'  

    script:
    """
    qiime taxa barplot \\
        --i-table ${feature_sample} \\
        --i-taxonomy ${taxonomy_qza} \\
        --m-metadata-file ${metadata} \\
        --o-visualization taxonomy_barplot.qzv
    """
}

