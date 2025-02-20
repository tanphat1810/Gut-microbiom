## Nextflow cho 16srRNA full-lenght tren gut microbiom

## Tổng quan

Pipeline Nextflow này xử lý dữ liệu microbiome sử dụng các định dạng đầu vào của QIIME 2. Nó thực hiện:

* Tóm tắt kiểm tra chất lượng (QC)

* Phân loại taxonomy bằng SILVA

* Trực quan hóa dữ liệu

Pipeline này được thiết kế cho các nghiên cứu phân tích microbiome, tập trung vào Gut microbiome

## Yêu cầu

* Nextflow

* Docker hoặc Singularity (để chạy container hóa)
* Java 8 hoặc cao hơn
## Cấu trúc thư mục mẫu

📦 **Gut-microbiome** (thư mục chứa các file làm việc, cần cd vào)
┣ 📜 main.nf
┣ 📜 nextflow.config
┣ 📜 workflow.nf
┣ 📂 conf
┃ ┗ 📜 base.config
┣ 📂 modules _(Lưu các module)_
┣ 📂 qiime_out (lưu đầu ra kết quả qiime2)
┣ 📂 fastqc_result (lưu đầu ra fastqc)
┗ 📂 seqkit_result (lưu đầu ra seqkit)
## Lệnh chạy nextflow

Lệnh chạy lúc bắt đầu 

nextflow run main.nf --input "data/*.fastq.gz" --outdir results -profile docker -c nextflow.config (không chỉ định params.json, đối với gộp params vào nextflow.config)

nextflow run main.nf --input "data/*.fastq.gz" --outdir results --params-file params.json --profile docker -c nextflow.config (đối với tách riêng params.json ra)


Lệnh resume

nextflow run main.nf --input "data/*.fastq.gz" --outdir results -profile docker -c nextflow.config -resume

