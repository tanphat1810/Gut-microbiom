#!/usr/bin/env nextflow
params.input = file("data/*.fastq.gz") // Đường dẫn file FASTQ
// Kiểm tra có file không
if (!params.input.size()) {
    println "ERROR: No FASTQ files found in the 'data/' directory. Please add files and try again."
    System.exit(1)
} else {
    println "FASTQ files detected: ${params.input}"
}
// Chạy FastQC
process run_fastqc {
    input:
    path fastq_files from params.input
    output:
    path "*.html"
    path "*.zip"
    script:
    """
    fastqc $fastq_files
    """
}
// Chạy MultiQC
process run_multiqc {
    input:
    path fastqc_zip.collect() 
    output:
    path "multiqc_report.html" 

    script:
    """
    multiqc .
    """
}
// cắt primer adapter lọc Q=20
params.input = file("data/*.fastq.gz")   // Đường dẫn dữ liệu FASTQ
params.adapters = file("adapters.fa")   // Đường dẫn trình tự adapter
params.minlen = 50                      // Độ dài tối thiểu sau lọc
params.qscore = 20                      // Ngưỡng Q-score tối thiểu
// Lọc Q-score bằng Trimmomatic
process trim_quality {
    input:
    path fastq_files from params.input.collect() 

    output:
    path "*.filtered.fastq.gz" into filtered_reads 

    script:
    """
    trimmomatic PE ${fastq_files[0]} ${fastq_files[1]} \
        R1_paired.filtered.fastq.gz R1_unpaired.filtered.fastq.gz \
        R2_paired.filtered.fastq.gz R2_unpaired.filtered.fastq.gz \
        SLIDINGWINDOW:4:${params.qscore} MINLEN:${params.minlen}
    """
}

// Cắt adapter bằng Cutadapt
process cut_adapters {
    input:
    path filtered_fastq from filtered_reads.collect() // Lấy dữ liệu sau lọc

    output:
    path "*.trimmed.fastq.gz" into trimmed_reads // Output FASTQ sau cắt adapter

    script:
    """
    cutadapt -a file:${params.adapters} -A file:${params.adapters} \
        -o R1_trimmed.fastq.gz -p R2_trimmed.fastq.gz \
        ${filtered_fastq[0]} ${filtered_fastq[1]}
    """
}
// Denoising into ASVs with DADA2
