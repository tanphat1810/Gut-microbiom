nextflow.enable.dsl=2

params.fastq_files = 'data/*.fastq'
params.silva_classifier = 'silva-138-99-classifier.qza'

process readsQC {
    container 'quay.io/biocontainers/fastqc:0.11.9--0'

    input:
    path fastq_file

    output:
    path 'fastqc_results/original/*', emit: original_qc

    script:
    """
    mkdir -p fastqc_results/original
    fastqc --outdir fastqc_results/original ${fastq_file}
    """
}

process multiqcOriginal {
    container 'quay.io/biocontainers/multiqc:1.12--pyhdfd78af_0'

    input:
    path 'fastqc_results/original/*'

    output:
    path 'fastqc_results/original/multiqc_report.html', emit: original_multiqc

    script:
    """
    multiqc -o fastqc_results/original fastqc_results/original
    """
}

process filterReadsFromData {
    container 'quay.io/biocontainers/seqkit:2.3.0--h9ee0642_0'

    input:
    path fastq_file

    output:
    path 'filtered_fastq/*', emit: filtered_reads

    script:
    """
    mkdir -p filtered_fastq
    seqkit seq -Q 20 -o filtered_fastq/filtered_${fastq_file.baseName}.fastq ${fastq_file}
    """
}

process readsQCFiltered {
    container 'quay.io/biocontainers/fastqc:0.11.9--0'

    input:
    path filtered_fastq_file

    output:
    path 'fastqc_results/filtered/*', emit: filtered_qc

    script:
    """
    mkdir -p fastqc_results/filtered
    fastqc --outdir fastqc_results/filtered ${filtered_fastq_file}
    """
}

process multiqcFiltered {
    container 'quay.io/biocontainers/multiqc:1.12--pyhdfd78af_0'

    input:
    path 'fastqc_results/filtered/*'

    output:
    path 'fastqc_results/filtered_multiqc_report.html', emit: filtered_multiqc

    script:
    """
    multiqc -o fastqc_results/filtered fastqc_results/filtered
    """
}

process importToQiime {
    container 'quay.io/qiime2/core:2023.2'

    input:
    path fastq_dir

    output:
    path 'imported_data/demux.qza', emit: demux_qza

    script:
    """
    mkdir -p imported_data
    qiime tools import \
        --type 'SampleData[PairedEndSequencesWithQuality]' \
        --input-path ${fastq_dir} \
        --input-format CasavaOneEightSingleLanePerSampleDirFmt \
        --output-path imported_data/demux.qza
    """
}

process denoiseCCS {
    container 'quay.io/qiime2/core:2023.2'

    input:
    path demux_qza

    output:
    path 'dada2_output/table.qza', emit: table_qza
    path 'dada2_output/rep_seqs.qza', emit: rep_seqs_qza
    path 'dada2_output/stats.qza', emit: stats_qza

    script:
    """
    mkdir -p dada2_output
    qiime dada2 denoise-ccs \
        --i-demultiplexed-seqs ${demux_qza} \
        --o-table dada2_output/table.qza \
        --o-representative-sequences dada2_output/rep_seqs.qza \
        --o-denoising-stats dada2_output/stats.qza
    """
}

process filterASVs {
    container 'quay.io/qiime2/core:2023.2'

    input:
    path table_qza

    output:
    path 'filtered_output/filtered_table.qza', emit: filtered_table_qza

    script:
    """
    mkdir -p filtered_output
    qiime feature-table filter-features \
        --i-table ${table_qza} \
        --p-min-frequency 5 \
        --o-filtered-table filtered_output/filtered_table.qza
    """
}

process removeChimeras {
    container 'quay.io/qiime2/core:2023.2'

    input:
    path filtered_table_qza

    output:
    path 'chimera_free_output/chimera_free_table.qza', emit: chimera_free_table_qza

    script:
    """
    mkdir -p chimera_free_output
    qiime vsearch uchime-denovo \
        --i-table ${filtered_table_qza} \
        --o-nonchimeras chimera_free_output/chimera_free_table.qza
    """
}

process assignTaxonomy {
    container 'quay.io/qiime2/core:2023.2'

    input:
    path rep_seqs_qza
    path params.silva_classifier

    output:
    path 'taxonomy_output/taxonomy.qza', emit: taxonomy_qza

    script:
    """
    mkdir -p taxonomy_output
    qiime feature-classifier classify-sklearn \
        --i-reads ${rep_seqs_qza} \
        --i-classifier ${params.silva_classifier} \
        --o-classification taxonomy_output/taxonomy.qza
    """
}

process generateBarplot {
    container 'quay.io/qiime2/core:2023.2'

    input:
    path taxonomy_qza
    path chimera_free_table_qza

    output:
    path 'barplot_output/taxonomy_barplot.qzv', emit: taxonomy_barplot_qzv

    script:
    """
    mkdir -p barplot_output
    qiime taxa barplot \
        --i-table ${chimera_free_table_qza} \
        --i-taxonomy ${taxonomy_qza} \
        --o-visualization barplot_output/taxonomy_barplot.qzv
    """
}

process buildPhyloTree {
    container 'quay.io/qiime2/core:2023.2'

    input:
    path rep_seqs_qza

    output:
    path 'phylo_output/rooted_tree.qza', emit: rooted_tree_qza

    script:
    """
    mkdir -p phylo_output
    qiime phylogeny align-to-tree-mafft-fasttree \
        --i-sequences ${rep_seqs_qza} \
        --o-rooted-tree phylo_output/rooted_tree.qza
    """
}

workflow {
    fastq_files = Channel.fromPath(params.fastq_files)

    qc_results = readsQC(fastq_files)
    multiqcOriginal(qc_results.original_qc)

    filtered_files = filterReadsFromData(fastq_files)

    qc_filtered = readsQCFiltered(filtered_files.filtered_reads)
    multiqcFiltered(qc_filtered.filtered_qc)

    // Gom tất cả các file FASTQ đã lọc vào một thư mục
    filtered_fastq_dir = filtered_files.filtered_reads.collectFile(name: 'filtered_fastq_dir')

    imported_data = importToQiime(filtered_fastq_dir)

    denoise_results = denoiseCCS(imported_data.demux_qza)
    filtered_asvs = filterASVs(denoise_results.table_qza)
    chimera_free_table = removeChimeras(filtered_asvs.filtered_table_qza)

    taxonomy = assignTaxonomy(denoise_results.rep_seqs_qza, params.silva_classifier)
    barplot = generateBarplot(taxonomy.taxonomy_qza, chimera_free_table.chimera_free_table_qza)

    buildPhyloTree(denoise_results.rep_seqs_qza)
}

