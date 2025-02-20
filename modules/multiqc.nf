process MULTIQC {
    label 'process_low'
    conda "bioconda::multiqc=1.14"
    container "quay.io/biocontainers/multiqc:1.14--pyhdfd78af_0"

    input:
    path (all_results)

    output:
    path "multiqc_output/multiqc_report.html", emit: html
    path "multiqc_output/multiqc_data", emit: data
   

    script:
    """
    mkdir -p multiqc_output
    multiqc --outdir multiqc_output .
    """
}

