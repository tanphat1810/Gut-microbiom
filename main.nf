#!/usr/bin/env nextflow
nextflow.enable.dsl=2
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PacBio HiFi 16S rRNA Nextflow Pipeline
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    - Kiểm tra chất lượng dữ liệu
    - Lọc dữ liệu
    - Denoising
    - Phân loại Taxonomy
    - Phân tích đa dạng sinh học
    - Trực quan hóa kết quả
----------------------------------------------------------------------------------------
*/
// **IMPORT WORKFLOW**
include { PACBIO_16S_PIPELINE } from './workflow.nf'
workflow {
    PACBIO_16S_PIPELINE()
}
