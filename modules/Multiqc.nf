process MULTIQC {
    label 'process_single'

    conda "bioconda::multiqc=1.14"
    container "quay.io/biocontainers/multiqc:1.14--pyhdfd78af_0"

    input:
    path fastqc_zips

    output:
    path "multiqc_report.html", emit: report
    path "multiqc_data", emit: data
    path "versions.yml", emit: versions

    script:
    """
    multiqc --title "FastQC Report" --outdir . $fastqc_zips
    echo 'MULTIQC_VERSION: ' \$(multiqc --version) > versions.yml
    """
}
