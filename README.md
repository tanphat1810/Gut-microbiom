# Nextflow cho 16srRNA full-lenght tren gut microbiom

# Tổng quan

Pipeline Nextflow này xử lý dữ liệu microbiome sử dụng các định dạng đầu vào của QIIME 2. Nó thực hiện:

* Tóm tắt kiểm tra chất lượng (QC)

* Phân loại taxonomy bằng SILVA

* Trực quan hóa dữ liệu

Pipeline này được thiết kế cho các nghiên cứu phân tích microbiome, tập trung vào Gut microbiome

# Yêu cầu

* Nextflow

* Docker hoặc Singularity (để chạy container hóa)
* Java 8 hoặc cao hơn
# Tệp đầu vào
## Tệp TSV mẫu
Tệp này phải bao gồm:
* sample-id: Mã định danh duy nhất cho mẫu
* absolute-filepath: Đường dẫn đầy đủ đến tệp dữ liệu cho từng mẫu
## Tệp TSV metadata
Tệp này phải bao gồm:
* sample_name: Tên của mẫu
* condition: Nhóm hoặc điều kiện mà mỗi mẫu thuộc về
## Sử dụng
* Tải pipeline.
* Chuẩn bị các tệp đầu vào (TSV mẫu và TSV metadata).
* Chạy pipeline bằng lệnh sau:
```bash
nextflow run main.nf --input samples.tsv --metadata metadata.tsv --outdir results
```
## Tham số
--input: Đường dẫn đến tệp TSV mẫu

--metadata: Đường dẫn đến tệp TSV metadata

--outdir: Thư mục nơi kết quả sẽ được lưu trữ

# Đầu ra

Pipeline tạo ra:

Tóm tắt QC: Các chỉ số kiểm tra chất lượng cho dữ liệu đầu vào

Kết quả phân loại Taxonomy: Phân loại SILVA cho từng mẫu

Trực quan hóa: Các biểu đồ đại diện cho kết quả
