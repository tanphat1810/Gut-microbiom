#!/usr/bin/env nextflow
nextflow.enable.dsl=2
def helpMessage() {
log.info"""
  Usage:
  This pipeline takes in the standard sample manifest and metadata file used in
  QIIME 2 and produces QC summary, taxonomy classification results and visualization.

  For samples TSV, two columns named "sample-id" and "absolute-filepath" are
  required. For metadata TSV file, at least two columns named "sample_name" and
  "condition" to separate samples into different groups.

  nextflow run main.nf --input samples.tsv --metadata metadata.tsv \\
    --dada2_cpu 8 --vsearch_cpu 8

  By default, sequences are first trimmed with cutadapt. If adapters are already trimmed, you can skip 
  cutadapt by specifying "--skip_primer_trim".

  Other important options:
  --front_p    Forward primer sequence. Default to F27. (default: AGRGTTYGATYMTGGCTCAG)
  --adapter_p    Reverse primer sequence. Default to R1492. (default: AAGTCGTAACAAGGTARCY)
  --filterQ    Filter input reads above this Q value (default: 20).
  --downsample    Limit reads to a maximum of N reads if there are more than N reads (default: off)
  --max_ee    DADA2 max_EE parameter. Reads with number of expected errors higher than
              this value will be discarded (default: 2)
  --minQ    DADA2 minQ parameter. Reads with any base lower than this score 
            will be removed (default: 0)
  --min_len    Minimum length of sequences to keep (default: 1000)
  --max_len    Maximum length of sequences to keep (default: 1600)
  --pooling_method    QIIME 2 pooling method for DADA2 denoise see QIIME 2 
                      documentation for more details (default: "pseudo", alternative: "independent") 
  --maxreject    max-reject parameter for VSEARCH taxonomy classification method in QIIME 2
                 (default: 100)
  --maxaccept    max-accept parameter for VSEARCH taxonomy classification method in QIIME 2
                 (default: 100)
  --min_asv_totalfreq    Total frequency of any ASV must be above this threshold
                         across all samples to be retained. Set this to 0 to disable filtering
                         (default 5)
  --min_asv_sample    ASV must exist in at least min_asv_sample to be retained. 
                      Set this to 0 to disable. (default 1)
  --vsearch_identity    Minimum identity to be considered as hit (default 0.97)
  --rarefaction_depth    Rarefaction curve "max-depth" parameter. By default the pipeline
                         automatically select a cut-off above the minimum of the denoised 
                         reads for >80% of the samples. This cut-off is stored in a file called
                         "rarefaction_depth_suggested.txt" file in the results folder
                         (default: null)
  --dada2_cpu    Number of threads for DADA2 denoising (default: 8)
  --vsearch_cpu    Number of threads for VSEARCH taxonomy classification (default: 8)
  --cutadapt_cpu    Number of threads for primer removal using cutadapt (default: 16)
  --outdir    Output directory name (default: "results")
  --vsearch_db	Location of VSEARCH database (e.g. silva-138-99-seqs.qza can be
                downloaded from QIIME database)
  --vsearch_tax    Location of VSEARCH database taxonomy (e.g. silva-138-99-tax.qza can be
                   downloaded from QIIME database)
  --silva_db   Location of Silva 138 database for taxonomy classification 
  --gtdb_db    Location of GTDB r202 for taxonomy classification
  --refseq_db    Location of RefSeq+RDP database for taxonomy classification
  --skip_primer_trim    Skip all primers trimming (switch off cutadapt and DADA2 primers
                        removal) (default: trim with cutadapt)
  --skip_nb    Skip Naive-Bayes classification (only uses VSEARCH) (default: false)
  --colorby    Columns in metadata TSV file to use for coloring the MDS plot
               in HTML report (default: condition)
  --run_picrust2    Run PICRUSt2 pipeline. Note that pathway inference with 16S using PICRUSt2
                    has not been tested systematically (default: false)
  --download_db    Download databases needed for taxonomy classification only. Will not
                   run the pipeline. Databases will be downloaded to a folder "databases"
                   in the Nextflow pipeline directory.
  --publish_dir_mode    Outputs mode based on Nextflow "publishDir" directive. Specify "copy"
                        if requires hard copies. (default: symlink)
  --version    Output version
  """
}
// Show help message
params.help = false
if (params.help) exit 0, helpMessage()
params.version = false
version = "0.7"
if (params.version) exit 0, log.info("$version")
params.download_db = false
params.skip_primer_trim = false
params.skip_nb = false
params.run_picrust2 = false
params.filterQ = 20
params.min_len = 1000
params.max_len = 1600
params.colorby = "condition"
params.skip_phylotree = false
params.minQ = 0

// Check input/
params.input = false
if (params.input){
  n_sample=file(params.input).countLines() - 1
  if (n_sample == 1) {
    params.min_asv_totalfreq = 0
    params.min_asv_sample = 0
    println("Only 1 sample. min_asv_sample and min_asv_totalfreq set to 0.")
  } else {
    params.min_asv_totalfreq = 5
    params.min_asv_sample = 1
  }
} else {
  println("No input file given to --input!")
  n_sample=0
  params.min_asv_totalfreq = 0
  params.min_asv_sample = 0
}

if (params.skip_primer_trim) {
  params.front_p = 'none'
  params.adapter_p = 'none'
  trim_cutadapt = "No"
} else {
  trim_cutadapt = "Yes"
// These are default V1-V9 adapter
params.front_p = 'AGRGTTYGATYMTGGCTCAG'
params.adapter_p = 'AAGTCGTAACAAGGTARCY'
process denoise_dada2 {
    tag "${sample_id}"

    input:
    path input_fastq, type: 'file' // FASTQ files cho từng mẫu
    path metadata, type: 'file'   // Tệp metadata TSV
    val sample_id                 // Tên mẫu (nếu cần)

    output:
    path "${sample_id}_ASVs.qza", emit: result // Tệp kết quả ASV QIIME2

    script:
    """
    qiime dada2 denoise-single \
        --i-demultiplexed-seqs ${input_fastq} \
        --p-trunc-len ${params.max_len} \
        --p-min-trunc-len ${params.min_len} \
        --p-max-ee ${params.max_ee} \
        --p-min-quality ${params.minQ} \
        --o-table ${sample_id}_ASVs_table.qza \
        --o-representative-sequences ${sample_id}_ASVs_rep_seqs.qza \
        --o-denoising-stats ${sample_id}_ASVs_denoising_stats.qza \
        --p-n-threads ${params.dada2_cpu}
    """
}

params {
    max_len = 250
    min_len = 100
    max_ee = 2
    minQ = 20
    dada2_cpu = 8
}
nextflow run main.nf --input input_manifest.tsv --metadata metadata.tsv
Channel.fromPath('input/*.fastq')
    .map { file -> tuple(file.baseName, file) }
    .set { fastq_files }

fastq_files
    .map { sample_id, input_fastq -> tuple(sample_id, input_fastq, metadata) }
    .into { dada2_input }

dada2_input
    .into(dada2_results)
    .process(dada2_results) // Có thể thêm bước xử lý tiếp theo.
Channel.fromPath('input/*.fastq')
    .map { file -> tuple(file.baseName, file) }
    .set { fastq_files }
nextflow.enable.dsl=2

params {
    input = null             // Tệp manifest đầu vào
    metadata = null          // Tệp metadata TSV
    max_len = 250            // Độ dài tối đa của trình tự sau khi cắt
    min_len = 100            // Độ dài tối thiểu của trình tự sau khi cắt
    max_ee = 2               // Số lỗi dự kiến tối đa
    minQ = 20                // Chất lượng tối thiểu
    dada2_cpu = 8            // Số luồng CPU cho DADA2
    outdir = 'results'       // Thư mục xuất kết quả
    min_asv_totalfreq = 5    // Tần suất tối thiểu cho ASVs
    min_asv_sample = 1       // Số lượng mẫu tối thiểu mà ASVs phải xuất hiện
}

workflow {
    // Bước đọc dữ liệu FASTQ từ thư mục input
    Channel.fromPath('input/*.fastq')
        .map { file -> tuple(file.baseName, file) }
        .set { fastq_files }

    // Truyền metadata vào từng mẫu FASTQ
    fastq_files
        .map { sample_id, input_fastq -> tuple(sample_id, input_fastq, file(params.metadata)) }
        .set { dada2_input }

    // Chạy bước DADA2 denoising
    process_denoise_dada2(dada2_input)

    // Bước lọc ASVs theo tần suất và mẫu
    dada2_results
        .set { filtered_asvs }
    process_filter_asvs(filtered_asvs)
}

// Process xử lý Denoising bằng DADA2
process process_denoise_dada2 {
    tag "${sample_id}"

    input:
    val sample_id                 // Tên mẫu
    path input_fastq, type: 'file' // Tệp FASTQ đầu vào cho từng mẫu
    path metadata, type: 'file'   // Metadata TSV đầu vào

    output:
    path "${sample_id}_ASVs_table.qza", emit: table // Bảng ASV
    path "${sample_id}_ASVs_rep_seqs.qza", emit: rep_seqs // Các trình tự đại diện
    path "${sample_id}_ASVs_denoising_stats.qza", emit: stats // Thống kê denoising

    script:
    """
    qiime dada2 denoise-single \
        --i-demultiplexed-seqs ${input_fastq} \
        --p-trunc-len ${params.max_len} \
        --p-min-length ${params.min_len} \
        --p-max-ee ${params.max_ee} \
        --p-min-quality ${params.minQ} \
        --o-table ${sample_id}_ASVs_table.qza \
        --o-representative-sequences ${sample_id}_ASVs_rep_seqs.qza \
        --o-denoising-stats ${sample_id}_ASVs_denoising_stats.qza \
        --p-n-threads ${params.dada2_cpu}
    """
}

// Process lọc ASVs theo tần suất tối thiểu và số mẫu tối thiểu
process process_filter_asvs {
    input:
    path table, type: 'file' // Bảng ASVs từ DADA2

    output:
    path "filtered_${table.baseName}", emit: filtered_table // Bảng ASVs sau khi lọc

    script:
    """
    qiime feature-table filter-features \
        --i-table ${table} \
        --p-min-frequency ${params.min_asv_totalfreq} \
        --p-min-samples ${params.min_asv_sample} \
        --o-filtered-table filtered_${table.baseName}
    """
}
nextflow run main.nf --input input_manifest.tsv --metadata metadata.tsv
nextflow.enable.dsl=2

params {
    input = null             // Tệp manifest đầu vào
    metadata = null          // Tệp metadata TSV
    max_len = 250            // Độ dài tối đa của trình tự sau khi cắt
    min_len = 100            // Độ dài tối thiểu của trình tự sau khi cắt
    max_ee = 2               // Số lỗi dự kiến tối đa
    minQ = 20                // Chất lượng tối thiểu
    dada2_cpu = 8            // Số luồng CPU cho DADA2
    outdir = 'results'       // Thư mục xuất kết quả
    min_asv_totalfreq = 5    // Tần suất tối thiểu cho ASVs
    min_asv_sample = 1       // Số lượng mẫu tối thiểu mà ASVs phải xuất hiện
}

workflow {
    // Bước đọc dữ liệu FASTQ từ thư mục input
    Channel.fromPath('input/*.fastq')
        .map { file -> tuple(file.baseName, file) }
        .set { fastq_files }

    // Truyền metadata vào từng mẫu FASTQ
    fastq_files
        .map { sample_id, input_fastq -> tuple(sample_id, input_fastq, file(params.metadata)) }
        .set { dada2_input }

    // Chạy bước DADA2 denoising
    process_denoise_dada2(dada2_input)

    // Bước lọc ASVs theo tần suất và mẫu
    dada2_results
        .set { filtered_asvs }
    process_filter_asvs(filtered_asvs)

    // Thu thập thống kê DADA2 từ kết quả denoising
    filtered_asvs
        .map { filtered_table -> tuple(filtered_table, file("${filtered_table.baseName}_ASVs_denoising_stats.qza")) }
        .set { dada2_stats }
    process_collect_dada2_stats(dada2_stats)
}

// Process xử lý Denoising bằng DADA2
process process_denoise_dada2 {
    tag "${sample_id}"

    input:
    val sample_id                 // Tên mẫu
    path input_fastq, type: 'file' // Tệp FASTQ đầu vào cho từng mẫu
    path metadata, type: 'file'   // Metadata TSV đầu vào

    output:
    path "${sample_id}_ASVs_table.qza", emit: table // Bảng ASV
    path "${sample_id}_ASVs_rep_seqs.qza", emit: rep_seqs // Các trình tự đại diện
    path "${sample_id}_ASVs_denoising_stats.qza", emit: stats // Thống kê denoising

    script:
    """
    qiime dada2 denoise-single \
        --i-demultiplexed-seqs ${input_fastq} \
        --p-trunc-len ${params.max_len} \
        --p-min-length ${params.min_len} \
        --p-max-ee ${params.max_ee} \
        --p-min-quality ${params.minQ} \
        --o-table ${sample_id}_ASVs_table.qza \
        --o-representative-sequences ${sample_id}_ASVs_rep_seqs.qza \
        --o-denoising-stats ${sample_id}_ASVs_denoising_stats.qza \
        --p-n-threads ${params.dada2_cpu}
    """
}

// Process lọc ASVs theo tần suất tối thiểu và số mẫu tối thiểu
process process_filter_asvs {
    input:
    path table, type: 'file' // Bảng ASVs từ DADA2

    output:
    path "filtered_${table.baseName}", emit: filtered_table // Bảng ASVs sau khi lọc

    script:
    """
    qiime feature-table filter-features \
        --i-table ${table} \
        --p-min-frequency ${params.min_asv_totalfreq} \
        --p-min-samples ${params.min_asv_sample} \
        --o-filtered-table filtered_${table.baseName}
    """
}

// Process thu thập thống kê DADA2
process process_collect_dada2_stats {
    input:
    path filtered_table, type: 'file' // Bảng ASVs đã lọc
    path stats, type: 'file'           // Tệp thống kê DADA2

    output:
    path "collected_denoising_stats.txt", emit: collected_stats

    script:
    """
    # Chuyển đổi và thu thập thống kê DADA2
    qiime metadata tabulate \
        --m-input-file ${stats} \
        --o-visualization denoising_stats_${filtered_table.baseName}.qzv

    # Lưu thông tin thống kê vào một tệp text
    qiime metadata tabulate \
        --m-input-file ${stats} \
        --o-visualization ${filtered_table.baseName}_denoising_stats.qzv

    # Tạo tệp text với thông tin thống kê
    echo "Denoising Stats for ${filtered_table.baseName}:" > collected_denoising_stats.txt
    qiime tools export --input-path denoising_stats_${filtered_table.baseName}.qzv --output-path denoising_stats_${filtered_table.baseName}
    cat denoising_stats_${filtered_table.baseName}/data/summary.tsv >> collected_denoising_stats.txt
    """
}
nextflow run main.nf --input input_manifest.tsv --metadata metadata.tsv

nextflow.enable.dsl=2

params {
    input = null             // Tệp manifest đầu vào
    metadata = null          // Tệp metadata TSV
    max_len = 250            // Độ dài tối đa của trình tự sau khi cắt
    min_len = 100            // Độ dài tối thiểu của trình tự sau khi cắt
    max_ee = 2               // Số lỗi dự kiến tối đa
    minQ = 20                // Chất lượng tối thiểu
    dada2_cpu = 8            // Số luồng CPU cho DADA2
    outdir = 'results'       // Thư mục xuất kết quả
    min_asv_totalfreq = 5    // Tần suất tối thiểu cho ASVs
    min_asv_sample = 1       // Số lượng mẫu tối thiểu mà ASVs phải xuất hiện
    vsearch_db = "GTDB_ssu_all_r207.qza"  // Cơ sở dữ liệu VSEARCH (có thể linh hoạt)
    vsearch_tax = "GTDB_ssu_all_r207.taxonomy.qza" // Tệp phân loại taxonomy
    skip_nb = false          // Có bỏ qua Naive-Bayes classifier không
}

workflow {
    // Đọc dữ liệu FASTQ từ thư mục input
    Channel.fromPath('input/*.fastq')
        .map { file -> tuple(file.baseName, file) }
        .set { fastq_files }

    // Truyền metadata vào từng mẫu FASTQ
    fastq_files
        .map { sample_id, input_fastq -> tuple(sample_id, input_fastq, file(params.metadata)) }
        .set { dada2_input }

    // Chạy bước DADA2 denoising
    process_denoise_dada2(dada2_input)

    // Lọc ASVs theo tần suất và mẫu
    dada2_results
        .set { filtered_asvs }
    process_filter_asvs(filtered_asvs)

    // Phân loại taxonomy với VSEARCH
    filtered_asvs
        .set { filtered_table }
    process_taxonomy_vsearch(filtered_table)

    // Phân loại taxonomy với Naive-Bayes (nếu không bỏ qua)
    if (!params.skip_nb) {
        filtered_asvs
            .set { filtered_table }
        process_taxonomy_nb(filtered_table)
    }
}

// Process xử lý Denoising bằng DADA2
process process_denoise_dada2 {
    tag "${sample_id}"

    input:
    val sample_id                 // Tên mẫu
    path input_fastq, type: 'file' // Tệp FASTQ đầu vào cho từng mẫu
    path metadata, type: 'file'   // Metadata TSV đầu vào

    output:
    path "${sample_id}_ASVs_table.qza", emit: table // Bảng ASV
    path "${sample_id}_ASVs_rep_seqs.qza", emit: rep_seqs // Các trình tự đại diện
    path "${sample_id}_ASVs_denoising_stats.qza", emit: stats // Thống kê denoising

    script:
    """
    qiime dada2 denoise-single \
        --i-demultiplexed-seqs ${input_fastq} \
        --p-trunc-len ${params.max_len} \
        --p-min-length ${params.min_len} \
        --p-max-ee ${params.max_ee} \
        --p-min-quality ${params.minQ} \
        --o-table ${sample_id}_ASVs_table.qza \
        --o-representative-sequences ${sample_id}_ASVs_rep_seqs.qza \
        --o-denoising-stats ${sample_id}_ASVs_denoising_stats.qza \
        --p-n-threads ${params.dada2_cpu}
    """
}

// Process lọc ASVs theo tần suất tối thiểu và số mẫu tối thiểu
process process_filter_asvs {
    input:
    path table, type: 'file' // Bảng ASVs từ DADA2

    output:
    path "filtered_${table.baseName}", emit: filtered_table // Bảng ASVs sau khi lọc

    script:
    """
    qiime feature-table filter-features \
        --i-table ${table} \
        --p-min-frequency ${params.min_asv_totalfreq} \
        --p-min-samples ${params.min_asv_sample} \
        --o-filtered-table filtered_${table.baseName}
    """
}

// Process phân loại taxonomy với VSEARCH
process process_taxonomy_vsearch {
    input:
    path filtered_table, type: 'file' // Bảng ASVs đã lọc
    path vsearch_db, type: 'file'      // Cơ sở dữ liệu VSEARCH
    path vsearch_tax, type: 'file'     // Tệp phân loại VSEARCH

    output:
    path "taxonomy_${filtered_table.baseName}", emit: taxonomy_result // Kết quả phân loại taxonomy

    script:
    """
    qiime vsearch classify-sklearn \
        --i-query ${filtered_table} \
        --i-reference-reads ${vsearch_db} \
        --i-reference-taxonomy ${vsearch_tax} \
        --o-classification taxonomy_${filtered_table.baseName}
    """
}

// Process phân loại taxonomy với Naive-Bayes (nếu không bỏ qua)
process process_taxonomy_nb {
    input:
    path filtered_table, type: 'file' // Bảng ASVs đã lọc

    output:
    path "taxonomy_nb_${filtered_table.baseName}", emit: nb_taxonomy_result // Kết quả phân loại taxonomy Naive-Bayes

    script:
    """
    qiime feature-classifier classify-sklearn \
        --i-classifier ${params.nb_classifier} \
        --i-reads ${filtered_table} \
        --o-classification taxonomy_nb_${filtered_table.baseName}
    """
}
nextflow run main.nf --input input_manifest.tsv --metadata metadata.tsv
nextflow.enable.dsl=2

params {
    input = null             // Tệp manifest đầu vào
    metadata = null          // Tệp metadata TSV
    max_len = 250            // Độ dài tối đa của trình tự sau khi cắt
    min_len = 100            // Độ dài tối thiểu của trình tự sau khi cắt
    max_ee = 2               // Số lỗi dự kiến tối đa
    minQ = 20                // Chất lượng tối thiểu
    dada2_cpu = 8            // Số luồng CPU cho DADA2
    outdir = 'results'       // Thư mục xuất kết quả
    min_asv_totalfreq = 5    // Tần suất tối thiểu cho ASVs
    min_asv_sample = 1       // Số lượng mẫu tối thiểu mà ASVs phải xuất hiện
    rarefaction_depth = null // Độ sâu rarefaction (nếu có)
}

workflow {
    // Đọc dữ liệu FASTQ từ thư mục input
    Channel.fromPath('input/*.fastq')
        .map { file -> tuple(file.baseName, file) }
        .set { fastq_files }

    // Truyền metadata vào từng mẫu FASTQ
    fastq_files
        .map { sample_id, input_fastq -> tuple(sample_id, input_fastq, file(params.metadata)) }
        .set { dada2_input }

    // Chạy bước DADA2 denoising
    process_denoise_dada2(dada2_input)

    // Lọc ASVs theo tần suất tối thiểu và số mẫu tối thiểu
    dada2_results
        .set { filtered_asvs }
    process_filter_asvs(filtered_asvs)

    // Thực hiện bước Rarefaction curve (nếu cần)
    filtered_asvs
        .set { rarefied_data }
    process_rarefaction_curve(rarefied_data)

    // Tiến hành các bước phân loại taxonomy với VSEARCH hoặc Naive-Bayes classifier nếu cần
    filtered_asvs
        .set { filtered_table }
    process_taxonomy_vsearch(filtered_table)

    if (!params.skip_nb) {
        filtered_asvs
            .set { filtered_table }
        process_taxonomy_nb(filtered_table)
    }
}

// Process Rarefaction curve
process process_rarefaction_curve {
    input:
    path filtered_table, type: 'file' // Bảng ASVs đã lọc

    output:
    path "rarefaction_curve_${filtered_table.baseName}.png", emit: rarefaction_plot // Đường cong rarefaction

    script:
    """
    qiime diversity alpha-rarefaction \
        --i-table ${filtered_table} \
        --p-max-depth ${params.rarefaction_depth} \
        --m-metadata-file ${params.metadata} \
        --o-visualization rarefaction_curve_${filtered_table.baseName}.qzv
    """
}

// Process Denoising bằng DADA2 (như trước)
process process_denoise_dada2 {
    tag "${sample_id}"

    input:
    val sample_id                 // Tên mẫu
    path input_fastq, type: 'file' // Tệp FASTQ đầu vào cho từng mẫu
    path metadata, type: 'file'   // Metadata TSV đầu vào

    output:
    path "${sample_id}_ASVs_table.qza", emit: table // Bảng ASV
    path "${sample_id}_ASVs_rep_seqs.qza", emit: rep_seqs // Các trình tự đại diện
    path "${sample_id}_ASVs_denoising_stats.qza", emit: stats // Thống kê denoising

    script:
    """
    qiime dada2 denoise-single \
        --i-demultiplexed-seqs ${input_fastq} \
        --p-trunc-len ${params.max_len} \
        --p-min-length ${params.min_len} \
        --p-max-ee ${params.max_ee} \
        --p-min-quality ${params.minQ} \
        --o-table ${sample_id}_ASVs_table.qza \
        --o-representative-sequences ${sample_id}_ASVs_rep_seqs.qza \
        --o-denoising-stats ${sample_id}_ASVs_denoising_stats.qza \
        --p-n-threads ${params.dada2_cpu}
    """
}

// Process lọc ASVs theo tần suất tối thiểu và số mẫu tối thiểu
process process_filter_asvs {
    input:
    path table, type: 'file' // Bảng ASVs từ DADA2

    output:
    path "filtered_${table.baseName}", emit: filtered_table // Bảng ASVs sau khi lọc

    script:
    """
    qiime feature-table filter-features \
        --i-table ${table} \
        --p-min-frequency ${params.min_asv_totalfreq} \
        --p-min-samples ${params.min_asv_sample} \
        --o-filtered-table filtered_${table.baseName}
    """
}

// Process phân loại taxonomy với VSEARCH
process process_taxonomy_vsearch {
    input:
    path filtered_table, type: 'file' // Bảng ASVs đã lọc
    path vsearch_db, type: 'file'      // Cơ sở dữ liệu VSEARCH
    path vsearch_tax, type: 'file'     // Tệp phân loại VSEARCH

    output:
    path "taxonomy_${filtered_table.baseName}", emit: taxonomy_result // Kết quả phân loại taxonomy

    script:
    """
    qiime vsearch classify-sklearn \
        --i-query ${filtered_table} \
        --i-reference-reads ${vsearch_db} \
        --i-reference-taxonomy ${vsearch_tax} \
        --o-classification taxonomy_${filtered_table.baseName}
    """
}

// Process phân loại taxonomy với Naive-Bayes (nếu không bỏ qua)
process process_taxonomy_nb {
    input:
    path filtered_table, type: 'file' // Bảng ASVs đã lọc

    output:
    path "taxonomy_nb_${filtered_table.baseName}", emit: nb_taxonomy_result // Kết quả phân loại taxonomy Naive-Bayes

    script:
    """
    qiime feature-classifier classify-sklearn \
        --i-classifier ${params.nb_classifier} \
        --i-reads ${filtered_table} \
        --o-classification taxonomy_nb_${filtered_table.baseName}
    """
}
// Process for Krona Visualization
process taxonomy_krona {
    input:
    path taxonomy_result, type: 'file' // Kết quả phân loại từ VSEARCH hoặc Naive-Bayes
    path metadata, type: 'file'        // Metadata TSV đầu vào

    output:
    path "krona_${taxonomy_result.baseName}.html", emit: krona_plot // Biểu đồ Krona

    script:
    """
    qiime taxa collapse \
        --i-table ${taxonomy_result} \
        --i-taxonomy ${taxonomy_result} \
        --p-level 6 \
        --o-collapsed-table collapsed_table.qza

    qiime tools export \
        --input-path collapsed_table.qza \
        --output-path collapsed_table

    biom convert \
        -i collapsed_table/feature-table.biom \
        -o collapsed_table/feature-table.tsv \
        --to-tsv

    python -c "import pandas as pd; df = pd.read_csv('collapsed_table/feature-table.tsv', sep='\\t', skiprows=1); df.to_csv('krona_input.tsv', sep='\\t', index=False)"

    ktImportText krona_input.tsv -o krona_${taxonomy_result.baseName}.html
    """
}

// Process for Taxonomy Barplot
process taxonomy_barplot {
    input:
    path taxonomy_result, type: 'file' // Kết quả phân loại từ VSEARCH hoặc Naive-Bayes
    path metadata, type: 'file'        // Metadata TSV đầu vào

    output:
    path "taxonomy_barplot_${taxonomy_result.baseName}.qzv", emit: barplot // Biểu đồ Barplot QIIME 2

    script:
    """
    qiime taxa barplot \
        --i-table ${taxonomy_result} \
        --i-taxonomy ${taxonomy_result} \
        --m-metadata-file ${metadata} \
        --o-visualization taxonomy_barplot_${taxonomy_result.baseName}.qzv
    """
}
workflow {
    // Các bước trước (denoising, filtering, rarefaction, taxonomy classification)
    filtered_asvs
        .set { taxonomy_results }
    process_taxonomy_vsearch(taxonomy_results)
    process_taxonomy_nb(taxonomy_results)

    // Bước Krona visualization
    taxonomy_results
        .map { taxonomy_result -> tuple(taxonomy_result, params.metadata) }
        .into { krona_input }
    krona_input
        .process(taxonomy_krona)

    // Bước Barplot visualization
    taxonomy_results
        .map { taxonomy_result -> tuple(taxonomy_result, params.metadata) }
        .into { barplot_input }
    barplot_input
        .process(taxonomy_barplot)
}
nextflow run main.nf --input input_manifest.tsv --metadata metadata.tsv
process generate_html_report {
    input:
    path taxonomy_result, type: 'file'       // Kết quả phân loại từ VSEARCH hoặc Naive-Bayes
    path barplot_result, type: 'file'       // Kết quả barplot QIIME2
    path krona_result, type: 'file'         // Biểu đồ Krona HTML
    path metadata, type: 'file'             // Tệp metadata TSV đầu vào
    path denoising_stats, type: 'file'      // Thống kê từ DADA2

    output:
    path "analysis_report.html", emit: report // Báo cáo HTML

    script:
    """
    Rscript -e "rmarkdown::render('${params.rmd_vis_biom_script}', output_file='analysis_report.html', params=list(
        taxonomy_result='${taxonomy_result}',
        barplot_result='${barplot_result}',
        krona_result='${krona_result}',
        metadata='${metadata}',
        denoising_stats='${denoising_stats}'
    ))"
    """
}
---
title: "16S rRNA Analysis Report"
output: html_document
params:
  taxonomy_result: "taxonomy.qza"
  barplot_result: "barplot.qzv"
  krona_result: "krona.html"
  metadata: "metadata.tsv"
  denoising_stats: "denoising-stats.qza"
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = FALSE, warning = FALSE, message = FALSE)
library(qiime2R)
library(ggplot2)
library(DT)
metadata <- read.csv(params$metadata, sep = "\t")
datatable(metadata, options = list(pageLength = 5))
taxonomy <- read_qza(params$taxonomy_result)
print(taxonomy)
print(params$barplot_result)
cat('<a href="', params$krona_result, '">Open Krona Visualization</a>')
denoising_stats <- read_qza(params$denoising_stats)
datatable(denoising_stats$data)

---

### **Tích hợp Workflow**

Thêm bước này vào workflow Nextflow để tạo báo cáo:

```groovy
workflow {
    // Các bước trước (denoising, filtering, taxonomy classification)
    taxonomy_results
        .set { taxonomy_results }

    barplot_results
        .set { barplot_results }

    krona_results
        .set { krona_results }

    taxonomy_results
        .map { taxonomy_result, barplot_result, krona_result, metadata, denoising_stats ->
            tuple(taxonomy_result, barplot_result, krona_result, metadata, denoising_stats)
        }
        .process(generate_html_report)
}
nextflow run main.nf --input input_manifest.tsv --metadata metadata.tsv
qiime tools export \
    --input-path feature-table.qza \
    --output-path exported_biom/
process export_to_biom {
    input:
    path feature_table, type: 'file' // Bảng ASV hoặc OTU từ bước taxonomy
    path taxonomy_result, type: 'file' // Kết quả phân loại taxonomy
    path metadata, type: 'file'        // Metadata TSV đầu vào

    output:
    path "exported_biom/feature-table.biom", emit: biom
    path "exported_biom/taxonomy.tsv", emit: taxonomy

    script:
    """
    # Xuất bảng tính năng (feature table) sang BIOM
    qiime tools export \
        --input-path ${feature_table} \
        --output-path exported_biom/

    # Chuyển taxonomy thành TSV cho tương thích với BIOM
    qiime tools export \
        --input-path ${taxonomy_result} \
        --output-path exported_biom/

    # Đổi tên để rõ ràng
    mv exported_biom/feature-table.biom exported_biom/feature-table.biom
    mv exported_biom/taxonomy.tsv exported_biom/taxonomy.tsv
    """
}
workflow {
    taxonomy_results
        .map { feature_table, taxonomy_result, metadata ->
            tuple(feature_table, taxonomy_result, metadata)
        }
        .process(export_to_biom)
}
biom add-metadata \
    -i feature-table.biom \
    -o feature-table-with-metadata.biom \
    --sample-metadata-file metadata.tsv
nextflow run main.nf --input input_manifest.tsv --metadata metadata.tsv
qiime phylogeny align-to-tree-mafft-fasttree \
    --i-sequences rep_seqs.qza \
    --o-alignment aligned_rep_seqs.qza \
    --o-masked-alignment masked_aligned_rep_seqs.qza \
    --o-tree unrooted_tree.qza \
    --o-rooted-tree rooted_tree.qza
process phylogenetic_tree {
    input:
    path rep_seqs, type: 'file'

    output:
    path "rooted_tree.qza", emit: rooted_tree
    path "unrooted_tree.qza", emit: unrooted_tree

    script:
    """
    qiime phylogeny align-to-tree-mafft-fasttree \
        --i-sequences ${rep_seqs} \
        --o-alignment aligned_rep_seqs.qza \
        --o-masked-alignment masked_aligned_rep_seqs.qza \
        --o-tree unrooted_tree.qza \
        --o-rooted-tree rooted_tree.qza
    """
}
qiime feature-table summarize \
    --i-table feature_table.qza \
    --o-visualization table_summary.qzv
qiime diversity core-metrics-phylogenetic \
    --i-phylogeny rooted_tree.qza \
    --i-table feature_table.qza \
    --p-sampling-depth 2000 \  # Thay giá trị 2000 bằng rarefaction depth
    --m-metadata-file metadata.tsv \
    --output-dir diversity_metrics
process diversity_metrics {
    input:
    path feature_table, type: 'file'
    path rooted_tree, type: 'file'
    path metadata, type: 'file'
    val rarefaction_depth // Giá trị user-defined hoặc auto-calculated

    output:
    path "diversity_metrics/", emit: diversity_results

    script:
    """
    qiime diversity core-metrics-phylogenetic \
        --i-phylogeny ${rooted_tree} \
        --i-table ${feature_table} \
        --p-sampling-depth ${rarefaction_depth} \
        --m-metadata-file ${metadata} \
        --output-dir diversity_metrics
    """
}
workflow {
    // Phylogenetic tree
    denoised_representative_seqs
        .set { rep_seqs }

    rep_seqs
        .process(phylogenetic_tree)
        .set { rooted_tree }

    // Diversity metrics
    feature_table
        .combine(rooted_tree, metadata)
        .map { table, tree, meta -> tuple(table, tree, meta, params.rarefaction_depth) }
        .process(diversity_metrics)
}
nextflow run main.nf --input input_manifest.tsv --metadata metadata.tsv --rarefaction_depth 2000
qiime metadata tabulate \
    --m-input-file diversity_metrics/shannon_vector.qza \
    --m-input-file diversity_metrics/observed_otus_vector.qza \
    --o-visualization diversity_metrics/alpha_diversity_summary.qzv

qiime metadata tabulate \
    --m-input-file taxonomy.qza \
    --o-visualization taxonomy_summary.qzv
qiime taxa barplot \
    --i-table feature_table.qza \
    --i-taxonomy taxonomy.qza \
    --m-metadata-file metadata.tsv \
    --o-visualization taxa_barplot.qzv
process hmlt_report {
    input:
    path alpha_metrics, type: 'file' // Kết quả diversity metrics (Shannon, Observed OTUs)
    path taxonomy, type: 'file'      // Kết quả taxonomy
    path metadata, type: 'file'      // Metadata mẫu

    output:
    path "hmlt_report/", emit: hmlt_results

    script:
    """
    mkdir -p hmlt_report

    qiime metadata tabulate \
        --m-input-file ${alpha_metrics}/shannon_vector.qza \
        --m-input-file ${alpha_metrics}/observed_otus_vector.qza \
        --o-visualization hmlt_report/alpha_diversity_summary.qzv

    qiime metadata tabulate \
        --m-input-file ${taxonomy} \
        --o-visualization hmlt_report/taxonomy_summary.qzv

    qiime taxa barplot \
        --i-table feature_table.qza \
        --i-taxonomy ${taxonomy} \
        --m-metadata-file ${metadata} \
        --o-visualization hmlt_report/taxa_barplot.qzv
    """
}
workflow {
    // Collect diversity metrics and taxonomy
    diversity_metrics
        .combine(taxonomy, metadata)
        .process(hmlt_report)
        .set { hmlt_results }

    // Generate final report
    hmlt_results
        .view { println "HMLT Report and Visualization generated in: $it" }
}
qiime tools export \
    --input-path feature_table.qza \
    --output-path exported_data/

qiime tools export \
    --input-path taxonomy.qza \
    --output-path exported_data/

