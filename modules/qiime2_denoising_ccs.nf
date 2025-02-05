process QIIME2_DENOISE_CCS {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::qiime2=2023.7"
    container "quay.io/biocontainers/qiime2-core:2023.7"

    input:
    tuple val(meta), path(demux_qza)

    output:
    tuple val(meta), path("*.table.qza"), emit: table_qza
    tuple val(meta), path("*.representative_sequences.qza"), emit: rep_seqs_qza
    tuple val(meta), path("*.denoising_stats.qza"), emit: stats_qza
    path "versions.yml", emit: versions

    script:
    """
    qiime dada2 denoise-ccs \\
        --i-demultiplexed-seqs $demux_qza \\
        --p-front '${params.denoise_ccs_front}' \\
        --p-adapter '${params.denoise_ccs_adapter}' \\
        --p-max-mismatch ${params.denoise_ccs_max_mismatch} \\
        --p-indels ${params.denoise_ccs_indels} \\
        --p-trunc-len ${params.denoise_ccs_trunc_len} \\
        --p-trim-left ${params.denoise_ccs_trim_left} \\
        --p-max-ee ${params.denoise_ccs_max_ee} \\
        --p-trunc-q ${params.denoise_ccs_trunc_q} \\
        --p-min-len ${params.denoise_ccs_min_len} \\
        --p-max-len ${params.denoise_ccs_max_len} \\
        --p-pooling-method ${params.denoise_ccs_pooling_method} \\
        --p-chimera-method ${params.denoise_ccs_chimera_method} \\
        --p-min-fold-parent-over-abundance ${params.denoise_ccs_min_fold_parent} \\
        --p-allow-one-off ${params.denoise_ccs_allow_one_off} \\
        --p-n-threads $task.cpus \\
        --p-n-reads-learn ${params.denoise_ccs_n_reads_learn} \\
        --p-hashed-feature-ids ${params.denoise_ccs_hashed_feature_ids} \\
        --p-retain-all-samples ${params.denoise_ccs_retain_all_samples} \\
        --o-table table.qza \\
        --o-representative-sequences rep-seqs.qza \\
        --o-denoising-stats stats.qza

    echo 'QIIME2_VERSION: ' \$(qiime --version) > versions.yml
    """
}
