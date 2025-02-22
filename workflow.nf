#!/usr/bin/env nextflow
nextflow.enable.dsl=2
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PacBio HiFi 16S rRNA Nextflow Workflow
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
// **IMPORT MODULES**
include { FASTQC_RAW } from './modules/fastqc_raw.nf'
include { FASTQC_FILTERED } from './modules/fastqc_filtered.nf'
include { SEQKIT_FILTER } from './modules/Seqkit.nf'
include { MULTIQC_RAW } from './modules/multiqc_raw.nf'
include { MULTIQC_FILTERED } from './modules/multiqc_filtered.nf'
include { GENERATE_TSV } from './modules/generate_tsv.nf'
include { QIIME_IMPORT } from './modules/import_qiime2.nf'
include { DENOISE_CCS } from './modules/denoise-ccs.nf'
include { FILTER_FEATURES } from './modules/qiime2_feature_table.nf'
include { FILTER_SAMPLES } from './modules/qiime2_feature_sample.nf'
include { FILTER_SEQS } from './modules/qiime2_filterseqs.nf'
include { CLASSIFY_TAXONOMY } from './modules/classify_taxonomy.nf'
include { BARPLOT } from './modules/taxonomy_barplot.nf'
include { BUILD_PHYLOGENETIC_TREE } from './modules/build_phylogenetictree.nf'
include { ALPHA_METRICS } from './modules/alpha_metrics.nf'
include { FAITH_PD } from './modules/alpha_faith_pd.nf'
include { BETA_DIVERSITY } from './modules/beta_diversity.nf'
include { BETA_PCOA } from './modules/beta_pcoa.nf'
include { BETA_EMPEROR } from './modules/beta_emperor.nf'
channel_reads = Channel.fromPath(params.input)
                       .map { [meta: [id: it.baseName], file: it] }
workflow PACBIO_16S_PIPELINE {
    // Chạy FastQC trên dữ liệu gốc
    FASTQC_RAW(channel_reads)
    // Chạy SEQKIT_FILTER để lọc dữ liệu
    SEQKIT_FILTER(channel_reads)
    // Nhận dữ liệu đã lọc từ SEQKIT_FILTER (đảm bảo chờ SEQKIT_FILTER hoàn tất)
    filtered_reads = SEQKIT_FILTER.out.filtered
        .map { meta, file -> tuple(meta, file) } // Đảm bảo giữ metadata đúng
        .view()
    // Chạy FastQC trên dữ liệu đã lọc
    FASTQC_FILTERED(filtered_reads)
    multiqc_raw_input = FASTQC_RAW.out.zip.map { meta, zip_file -> zip_file }.collect()
    multiqc_raw_result = MULTIQC_RAW(multiqc_raw_input)
    multiqc_filtered_input = FASTQC_FILTERED.out.zip.map { meta, zip_file -> zip_file }.collect()
    multiqc_filtered_result = MULTIQC_FILTERED(multiqc_filtered_input)
filtered_reads_path = SEQKIT_FILTER.out.filtered
    .map { it[1] }  // Chỉ lấy đường dẫn file
    .collect()
    .view()
GENERATE_TSV(filtered_reads_path)
manifest = GENERATE_TSV.out
    .filter { it.name.endsWith('.tsv') }  // Chỉ lấy file có đuôi .tsv
    .flatten()
    .view()
QIIME_IMPORT(manifest)
    qiime_data = QIIME_IMPORT.out.qiime_artifact
    .toList()
    .view()
    DENOISE_CCS(qiime_data) 
     FILTER_FEATURES(DENOISE_CCS.out.feature_table)
     FILTER_SAMPLES(FILTER_FEATURES.out.filtered_feature_table)
     filtered_table = FILTER_FEATURES.out.filtered_feature_table
     FILTER_SEQS(
        DENOISE_CCS.out.representative_sequences,
        FILTER_SAMPLES.out.filtered_sample_table
    )
classifier_path = file(params.classifier)
CLASSIFY_TAXONOMY(
        FILTER_SEQS.out.filtered_representative_seqs,
        classifier_path
        )
metadata_file = file(params.metadata)  
BARPLOT(
        CLASSIFY_TAXONOMY.out.classified_taxonomy,
        FILTER_SAMPLES.out.filtered_sample_table,
        metadata_file
)
    BUILD_PHYLOGENETIC_TREE(FILTER_SEQS.out.filtered_representative_seqs)
    ALPHA_METRICS(FILTER_SAMPLES.out.filtered_sample_table)
    FAITH_PD(
        FILTER_SAMPLES.out.filtered_sample_table,
        BUILD_PHYLOGENETIC_TREE.out.rooted_tree)
    BETA_DIVERSITY(FILTER_SAMPLES.out.filtered_sample_table, BUILD_PHYLOGENETIC_TREE.out.rooted_tree) 
   BETA_PCOA(
    BETA_DIVERSITY.out.jaccard, 
    BETA_DIVERSITY.out.bray_curtis, 
    BETA_DIVERSITY.out.unweighted_unifrac, 
    BETA_DIVERSITY.out.weighted_unifrac
)
    BETA_EMPEROR(BETA_PCOA.out.pcoa_results, metadata_file)
}
