process DENOISE_CCS {
    tag "qiime2_denoise_ccs"

    conda "bioconda::qiime2=2023.7"
    container "quay.io/qiime2/amplicon:2024.5"

    input:
    path demux_qza  // Nhận tệp đầu vào từ QIIME_IMPORT

    output:
    path "table.qza", emit: feature_table
    path "rep_seqs.qza", emit: representative_sequences
    path "denoising_stats.qza", emit: denoising_stats
    path "table.qzv", emit: table_qzv
    path "rep_seqs.qzv", emit: rep_seqs_qzv
    path "denoising_stats.qzv", emit: denoising_stats_qzv
    publishDir "qiime_out", mode: 'copy'
 script:
    """
     qiime dada2 denoise-ccs \
   --i-demultiplexed-seqs $demux_qza \
   --p-front "${params.qiime_denoise.front}" \
   --p-adapter "${params.qiime_denoise.adapter}" \
   --p-max-mismatch ${params.qiime_denoise.max_mismatch} \
   --p-indels ${params.qiime_denoise.indels} \
   --p-trunc-len ${params.qiime_denoise.trunc_len} \
   --p-trim-left ${params.qiime_denoise.trim_left} \
   --p-max-ee ${params.qiime_denoise.max_ee} \
   --p-trunc-q ${params.qiime_denoise.trunc_q} \
   --p-min-len ${params.qiime_denoise.min_len} \
   --p-max-len ${params.qiime_denoise.max_len} \
   --p-pooling-method ${params.qiime_denoise.pooling_method} \
   --p-chimera-method ${params.qiime_denoise.chimera_method} \
   --p-min-fold-parent-over-abundance ${params.qiime_denoise.min_fold_parent_over_abundance} \
   --p-allow-one-off ${params.qiime_denoise.allow_one_off} \
   --p-n-threads ${params.qiime_denoise.n_threads} \
   --p-n-reads-learn ${params.qiime_denoise.n_reads_learn} \
   --p-hashed-feature-ids ${params.qiime_denoise.hashed_feature_ids} \
   --o-table table.qza \
   --o-representative-sequences rep_seqs.qza \
   --o-denoising-stats denoising_stats.qza
   qiime feature-table summarize \
    --i-table table.qza \
    --o-visualization table.qzv
   qiime feature-table tabulate-seqs \
    --i-data rep_seqs.qza \
    --o-visualization rep_seqs.qzv
   qiime metadata tabulate \
    --m-input-file denoising_stats.qza \
    --o-visualization denoising_stats.qzv
    """
}
