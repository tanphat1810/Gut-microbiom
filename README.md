## Nextflow cho 16srRNA full-lenght tren gut microbiom

## Tổng quan

Pipeline Nextflow này xử lý dữ liệu microbiome sử dụng các định dạng đầu vào của QIIME 2. Nó thực hiện:

Kiểm tra chất lượng đọc (Reads QC)

Dùng FastQC hoặc MultiQC để kiểm tra chất lượng các đoạn đọc.
Lọc dữ liệu các đoạn đọc > Q20 (Reads Filtering)

Sử dụng Seqkit để lọc các đoạn đọc có chất lượng thấp.
Loại bỏ nhiễu và cắt tỉa, tạo ASVs (Trim and Denoising)

Dùng QIIME 2 để loại bỏ nhiễu và tạo ASVs.
Lọc bỏ tần suất tối thiểu và mẫu tối thiểu

Thiết lập ngưỡng lọc cho các ASVs xuất hiện ít hơn số lần nhất định.
Lọc các ASVs nhiễm chéo (Chimeric ASVs Filtering)

Dùng QIIME 2 để loại bỏ các ASVs bị nhiễm chéo.
Phân loại taxon bằng SILVA Database

Sử dụng QIIME 2 để gán danh tính loài cho ASVs với độ tương đồng 97%.
Phân tích cây phát sinh loài và các chỉ số đa dạng

Dùng QIIME 2 để xây dựng cây phát sinh loài và tính các chỉ số đa dạng (alpha, beta diversity).
Tạo biểu đồ phân loại vi sinh vật

Dùng QIIME 2 view để tạo biểu đồ phân loại (taxonomy barplot).
Xuất báo cáo và trực quan hóa dữ liệu

Dùng R script để tạo báo cáo và biểu đồ trực quan.

## cấu trúc thư mục mẫu

📂 Gut-microbiom

┣ 📜 main.nf

┣ 📜 workflow.nf

┣ 📜 nextflow.config

┣ 📜 params.json

┣ 📂 conf

┃ ┗ 📜 base.config

┣ 📂 modules _(Lưu các module)_

┣ 📂 metadata (lưu file metadata cho tạo taxonomy, tạo trước khi chạy nextflow)

┣ 📂 qiime_out (kết quả các process sẽ được copy vào, tạo trước khi chạy nextflow)

┣ 📂 result (nextflow tự tạo, lưu tiến trình chạy nextflow)

┣ 📂 work (lưu kết quả chính của các process, nextflow tự tạo, có thể dùng để resume)

## Lệnh chạy nextflow

Lệnh chạy lúc bắt đầu 

nextflow run main.nf  --outdir results -profile docker -c nextflow.config 

(không chỉ định params.json, đối với gộp params vào nextflow.config)

nextflow run main.nf -params-file params.json -profile docker -c nextflow.config 

(đối với tách riêng params.json ra)


Lệnh resume

nextflow run main.nf --outdir results -profile docker -c nextflow.config -resume

