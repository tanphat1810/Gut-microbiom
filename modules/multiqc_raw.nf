process MULTIQC_RAW {
    label 'process_low'
    conda "bioconda::multiqc=1.14"
    container "quay.io/biocontainers/multiqc:1.14--pyhdfd78af_0"

    input:
    path (all_results)

    output:
    path "multiqc_output/multiqc_raw_report.html", emit: html
    path "multiqc_output/multiqc_raw_data", emit: data
    publishDir "qiime_out", mode: 'copy'

    script:
    """
    mkdir -p multiqc_output
    multiqc --filename multiqc_raw_report.html --outdir multiqc_output .
    mv multiqc_output/multiqc_raw_report_data multiqc_output/multiqc_raw_data
    """
}

