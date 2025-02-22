process QIIME_IMPORT {
    tag "qiime2_import"
    
    conda "bioconda::qiime2=2023.7"
    container "quay.io/qiime2/amplicon:2024.5"

    input:
    path manifest  // Nhận tệp TSV đầu vào

    output:
    path "demux.qza", emit: qiime_artifact  // Xuất ra tệp QZA của QIIME2
    path "demux_summary.qzv", emit: demux_qzv // Xuất ra visualization file
    publishDir "qiime_out", mode: 'copy'

    script:
    """
    qiime tools import \\
        --type 'SampleData[SequencesWithQuality]' \\
        --input-path ${manifest} \\
        --output-path demux.qza \\
        --input-format SingleEndFastqManifestPhred33V2
    qiime demux summarize \\
        --i-data demux.qza \\
        --o-visualization demux_summary.qzv
    """
}

