process CREATE_METADATA {
    tag "$meta.id"
    label 'process_single'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("metadata.tsv"), emit: metadata

    script:
    """
    echo -e "#SampleID\\tabsolute-filepath" > metadata.tsv
    for file in $reads; do
        sample_id=\$(basename \$file | sed 's/.fastq.gz//')
        echo -e "\$sample_id\\t\$PWD/\$file" >> metadata.tsv
    done
    """
}
