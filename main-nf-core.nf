#!/usr/bin/env nextflow

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Pipeline xử lý dữ liệu 16S rRNA từ PacBio HiFi
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Mô tả  : Xử lý reads 16S full-length từ PacBio HiFi
    Công cụ: FastQC, SeqKit, QIIME2, R
    Github : (Thêm link nếu có)
----------------------------------------------------------------------------------------
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { FASTQC                  } from './modules/nf-core/fastqc/main'
include { MULTIQC                 } from './modules/nf-core/multiqc/main'
include { SEQKIT_FILTER            } from './modules/local/seqkit_filter'
include { QIIME2_DENOISE           } from './modules/local/qiime2_denoise'
include { QIIME2_FILTER_MIN_FREQ   } from './modules/local/qiime2_filter_min_freq'
include { QIIME2_FILTER_CHIMERA    } from './modules/local/qiime2_filter_chimera'
include { QIIME2_TAXONOMY_SILVA    } from './modules/local/qiime2_taxonomy_silva'
include { QIIME2_TREE_DIVERSITY    } from './modules/local/qiime2_tree_diversity'
include { QIIME2_TAXONOMY_BAR_PLOT } from './modules/local/qiime2_taxonomy_barplot'
include { REPORT_VISUALIZATION     } from './modules/local/report_visualization'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NAMED WORKFLOWS FOR PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PACBIO_16S_PIPELINE {

    main:

    //
    // Kiểm tra chất lượng reads với FastQC + MultiQC
    //
    FASTQC ( params.input_reads )
    MULTIQC ( FASTQC.out.zip )

    //
    // Lọc reads có chất lượng Q > 20 với SeqKit
    //
    SEQKIT_FILTER ( FASTQC.out.fastq, params.q_min )

    //
    // Loại bỏ nhiễu và tạo ASVs với QIIME2 DADA2
    //
    QIIME2_DENOISE ( SEQKIT_FILTER.out.filtered_reads, params.trunc_len )

    //
    // Lọc ASVs theo tần suất tối thiểu và số mẫu tối thiểu
    //
    QIIME2_FILTER_MIN_FREQ ( QIIME2_DENOISE.out.asv_table, params.min_freq, params.min_samples )

    //
    // Lọc ASVs nhiễm chéo (chimera)
    //
    QIIME2_FILTER_CHIMERA ( QIIME2_FILTER_MIN_FREQ.out.filtered_asv_table )

    //
    // Phân loại taxon với CSDL SILVA (độ tương đồng 97%)
    //
    QIIME2_TAXONOMY_SILVA ( QIIME2_FILTER_CHIMERA.out.filtered_asv_table, params.silva_db )

    //
    // Xây dựng cây phát sinh loài & chỉ số đa dạng
    //
    QIIME2_TREE_DIVERSITY ( QIIME2_FILTER_CHIMERA.out.filtered_asv_table, params.metadata )

    //
    // Biểu đồ phân loại sinh vật (Taxonomy barplot)
    //
    QIIME2_TAXONOMY_BAR_PLOT ( QIIME2_TAXONOMY_SILVA.out.taxonomy_table )

    //
    // Xuất báo cáo & biểu đồ trực quan hoá
    //
    REPORT_VISUALIZATION (
        MULTIQC.out.report,
        QIIME2_TAXONOMY_BAR_PLOT.out.barplot,
        QIIME2_TREE_DIVERSITY.out.diversity_metrics
    )

    emit:
    multiqc_report = MULTIQC.out.report
}
