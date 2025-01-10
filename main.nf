
nextflow.enable.dsl=2

process readsQC {
    // FastQC và MultiQC để kiểm tra chất lượng đoạn đọc
    container 'biocontainers/fastqc:v0.11.9_cv8'
    container 'biocontainers/multiqc:v1.11_cv1'

    input:
    path fastq_files from params.fastq_files

    output:
    path "fastqc_results/original/*"
    path "fastqc_results/original_multiqc_report.html"

    script:
    """
    mkdir -p fastqc_results/original
    fastqc --outdir fastqc_results/original ${fastq_files}
    multiqc -o fastqc_results/original fastqc_results/original
    """
}

process filterReads {
    // SeqKit để lọc dữ liệu đoạn đọc Q=20
    container 'quay.io/biocontainers/seqkit:2.3.0--h9a82719_0'

    input:
    path fastq_files from readsQC.out.collect()

    output:
    path "filtered_fastq/*"

    script:
    """
    if [ ! -f ${fastq_files} ]; then
        echo "Error: Filtered FASTQ files not found!"
        exit 1
    fi

    mkdir -p filtered_fastq
    seqkit seq -q 20 -o filtered_fastq/filtered_${fastq_files} ${fastq_files}
    """
}

process readsQCFiltered {
    // FastQC và MultiQC để kiểm tra chất lượng đoạn đọc sau khi lọc
    container 'biocontainers/fastqc:v0.11.9_cv8'
    container 'biocontainers/multiqc:v1.11_cv1'

    input:
    path filtered_fastq_files from filterReads.out.collect()

    output:
    path "fastqc_results/filtered/*"
    path "fastqc_results/filtered_multiqc_report.html"

    script:
    """
    mkdir -p fastqc_results/filtered
    fastqc --outdir fastqc_results/filtered ${filtered_fastq_files}
    multiqc -o fastqc_results/filtered fastqc_results/filtered
    """
}

process importToQiime {
    // QIIME 2 để import dữ liệu
    container 'qiime2/core:2023.2'

    input:
    path fastq_files from filterReads.out.collect()

    output:
    path "imported_data/demux.qza"

    script:
    """
    if [ ! -f ${fastq_files} ]; then
        echo "Error: Filtered FASTQ files not found for import!"
        exit 1
    fi

    mkdir -p imported_data

    qiime tools import \
        --type 'SampleData[PairedEndSequencesWithQuality]' \
        --input-path ${fastq_files} \
        --input-format CasavaOneEightSingleLanePerSampleDirFmt \
        --output-path imported_data/demux.qza
    """
}

process denoiseCCS {
    // QIIME 2 để loại bỏ nhiễu (dùng denoise ccs)
    container 'qiime2/core:2023.2'

    input:
    path demux_qza from importToQiime.out

    output:
    path "dada2_output/table.qza"
    path "dada2_output/rep_seqs.qza"
    path "dada2_output/stats.qza"

    script:
    """
    if [ ! -f ${demux_qza} ]; then
        echo "Error: Demultiplexed sequences not found!"
        exit 1
    fi

    mkdir -p dada2_output

    qiime dada2 denoise-ccs \
        --i-demultiplexed-seqs ${demux_qza} \
        --o-table dada2_output/table.qza \
        --o-representative-sequences dada2_output/rep_seqs.qza \
        --o-denoising-stats dada2_output/stats.qza
    """
}

process filterASVs {
    // Lọc bỏ các ASVs và mẫu không đạt tần suất tối thiểu
    container 'qiime2/core:2023.2'

    input:
    path table_qza from denoiseCCS.out.collect()[0]

    output:
    path "filtered_output/filtered_table.qza"

    script:
    """
    if [ ! -f ${table_qza} ]; then
        echo "Error: ASV table not found!"
        exit 1
    fi

    mkdir -p filtered_output

    qiime feature-table filter-features \
        --i-table ${table_qza} \
        --p-min-frequency 5 \
        --o-filtered-table filtered_output/filtered_table.qza

    qiime feature-table filter-samples \
        --i-table filtered_output/filtered_table.qza \
        --p-min-features 1 \
        --o-filtered-table filtered_output/filtered_table.qza
    """
}

process removeChimeras {
    // loại bỏ các ASVs nhiễm chéo
    container 'qiime2/core:2023.2'

    input:
    path filtered_table from filterASVs.out.collect()

    output:
    path "chimera_free_output/chimera_free_table.qza"

    script:
    """
    if [ ! -f ${filtered_table} ]; then
        echo "Error: Filtered ASV table not found!"
        exit 1
    fi

    mkdir -p chimera_free_output

    qiime vsearch uchime-denovo \
        --i-table ${filtered_table} \
        --o-chimeras chimera_free_output/chimeras.qza \
        --o-nonchimeras chimera_free_output/chimera_free_table.qza
    """
}

process assignTaxonomy {
    // Cơ sở dữ liệu SILVA để phân loại taxon
    container 'qiime2/core:2023.2'

    input:
    path rep_seqs from denoiseCCS.out.collect()[1]

    output:
    path "taxonomy_output/taxonomy.qza"

    script:
    """
    if [ ! -f ${rep_seqs} ]; then
        echo "Error: Representative sequences not found!"
        exit 1
    fi

    mkdir -p taxonomy_output

    qiime feature-classifier classify-sklearn \
        --i-reads ${rep_seqs} \
        --i-classifier silva-138-99-klassifier.qza \
        --o-classification taxonomy_output/taxonomy.qza
    """
}

process generateBarplot {
    // tạo taxa plot
    container 'qiime2/core:2023.2'

    input:
    path taxonomy_qza from assignTaxonomy.out
    path filtered_table from removeChimeras.out.collect()

    output:
    path "barplot_output/taxonomy_barplot.qzv"

    script:
    """
    if [ ! -f ${taxonomy_qza} ] || [ ! -f ${filtered_table} ]; then
        echo "Error: Taxonomy file or filtered table not found!"
        exit 1
    fi

    mkdir -p barplot_output

    qiime taxa barplot \
        --i-table ${filtered_table} \
        --i-taxonomy ${taxonomy_qza} \
        --o-visualization barplot_output/taxonomy_barplot.qzv
    """
}

process buildPhyloTree {
    // xây dựng cây phát sinh loài và tính toán các chỉ số đa dạng
    container 'qiime2/core:2023.2'
    input:
    path rep_seqs from denoiseCCS.out.collect()[1]
    output:
    path "phylo_output/rooted_tree.qza"
    script:
    """
    if [ ! -f ${rep_seqs} ]; then
        echo "Error: Representative sequences not found for phylogenetic tree!"
        exit 1
    fi
    mkdir -p phylo_output

    qiime phylogeny align-to-tree-mafft-fasttree \
        --i-sequences ${rep_seqs} \
        --o-rooted-tree phylo_output/rooted_tree.qza
    """
}
workflow {
    Channel.fromPath(params.fastq_files).set { fastq_files }
    readsQC(fastq_files)
    filterReads(fastq_files)
    readsQCFiltered()
    importToQiime()
    denoiseCCS()
    filterASVs()
    removeChimeras()
