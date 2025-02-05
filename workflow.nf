#!/usr/bin/env nextflow

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PacBio HiFi 16S rRNA Nextflow Pipeline
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    - Kiểm tra chất lượng dữ liệu
    - Lọc dữ liệu
    - Nhóm ASVs
    - Loại bỏ nhiễu
    - Phân loại Taxonomy
    - Phân tích đa dạng sinh học
    - Trực quan hóa kết quả
----------------------------------------------------------------------------------------
*/

include { FASTQC_MULTIQC }            from './modules/local/fastqc_multiqc.nf'
include { SEQKIT_FILTER }              from './modules/local/seqkit_filter.nf'
include { CREATE_TSV }                 from './modules/local/create_tsv.nf'
include { IMPORT_QIIME2 }              from './modules/local/import_qiime2.nf'
include { DEMUX_VISUALIZE }            from './modules/local/qiime2_demux_visualize.nf'
include { DENOISE_CCS }                from './modules/local/qiime2_denoise_ccs.nf'
include { DADA2_SUMMARY }              from './modules/local/dada2_summary.nf'
include { FEATURE_TABLE_SUMMARY }      from './modules/local/feature_table_summary.nf'
include { FEATURE_TABLE_TABULATE }     from './modules/local/feature_table_tabulate.nf'
include { FILTER_FEATURES }            from './modules/local/filter_features.nf'
include { CHIMERA_FILTERING }          from './modules/local/chimera_filtering.nf'
include { TAXONOMY_CLASSIFICATION }    from './modules/local/taxonomy_classification.nf'
include { TAXONOMY_VISUALIZATION }     from './modules/local/taxonomy_visualization.nf'
include { PHYLOGENY_TREE }             from './modules/local/phylogeny_tree.nf'
include { ALPHA_DIVERSITY }            from './modules/local/alpha_diversity.nf'
include { BETA_DIVERSITY }             from './modules/local/beta_diversity.nf'
include { ALPHA_DIVERSITY_VISUALIZATION } from './modules/local/alpha_diversity_visualization.nf'
include { BETA_DIVERSITY_VISUALIZATION }  from './modules/local/beta_diversity_visualization.nf'

workflow PACBIO_16S_PIPELINE {

    //  **Bước 1: Kiểm tra chất lượng dữ liệu**
    Channel.fromPath(params.input)
        | FASTQC_MULTIQC()

    //  **Bước 2: Lọc reads Q20**
    FASTQC_MULTIQC.out.filtered_reads
        | SEQKIT_FILTER()

    //  **Bước 3: Tạo TSV & Import dữ liệu vào QIIME2**
    SEQKIT_FILTER.out.filtered_reads
        | CREATE_TSV()
        | IMPORT_QIIME2()

    //  **Bước 4: Kiểm tra demultiplexing**
    IMPORT_QIIME2.out.qza
        | DEMUX_VISUALIZE()

    //  **Bước 5: Denoising ASVs bằng DADA2 CCS**
    IMPORT_QIIME2.out.qza
        | DENOISE_CCS()

    //  **Bước 6: Kiểm tra kết quả DADA2**
    DENOISE_CCS.out.feature_table
        | DADA2_SUMMARY()
        | FEATURE_TABLE_SUMMARY()
        | FEATURE_TABLE_TABULATE()

    //  **Bước 7: Lọc ASVs tần suất thấp**
    DADA2_SUMMARY.out.feature_table
        | FILTER_FEATURES()

    //  **Bước 8: Loại bỏ ASVs nhiễm chéo (Chimera)**
    FILTER_FEATURES.out.filtered_table
        | CHIMERA_FILTERING()

    //  **Bước 9: Phân loại Taxonomy**
    CHIMERA_FILTERING.out.cleaned_table
        | TAXONOMY_CLASSIFICATION()

    //  **Bước 10: Báo cáo phân loại**
    TAXONOMY_CLASSIFICATION.out.taxonomy
        | TAXONOMY_VISUALIZATION()

    //  **Bước 11: Xây dựng cây phát sinh loài**
    CHIMERA_FILTERING.out.cleaned_table
        | PHYLOGENY_TREE()

    //  **Bước 12: Tính chỉ số đa dạng Alpha & Beta**
    CHIMERA_FILTERING.out.cleaned_table
        | ALPHA_DIVERSITY()
        | BETA_DIVERSITY()

    //  **Bước 13: Trực quan hóa chỉ số đa dạng**
    ALPHA_DIVERSITY.out.alpha_metrics
        | ALPHA_DIVERSITY_VISUALIZATION()
    BETA_DIVERSITY.out.beta_metrics
        | BETA_DIVERSITY_VISUALIZATION()
}
