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
- Gut-microbiom (tên có thể đặt tùy ý, đây là thư mục CD vào để chạy nextflow)
 -- main.nf
 -- nextflow.config
 -- workflow.nf
 -- params.json
 -- 

