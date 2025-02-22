process GENERATE_TSV {
    tag "Generate TSV for QIIME2"
    
    input:
    path filtered_reads

    output:
    path "qiime2_input.tsv"
    publishDir "qiime_out", mode: 'copy'
    script:
    """
    echo -e "sample-id\\tabsolute-filepath" > qiime2_input.tsv
    for file in ${filtered_reads}; do
        sample_id=\$(basename "\$file" | sed 's/.filtered.fastq.gz//')
        echo -e "\$sample_id\\t\$PWD/\$file" >> qiime2_input.tsv
    done
    """
}

