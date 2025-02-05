process IQTREE_PHYLOGENY {
    tag "$meta.id"
    label 'process_high'

    conda "bioconda::iqtree=2.2.2.6 bioconda::qiime2=2023.7"
    container "quay.io/biocontainers/iqtree:2.2.2.6--h56fc30b_0"

    input:
    tuple val(meta), path(aligned_rep_seqs)
    tuple val(meta), path(feature_table)

    output:
    tuple val(meta), path("rooted-tree.qza"), emit: rooted_tree
    tuple val(meta), path("unrooted-tree.qza"), emit: unrooted_tree
    path "rarefaction_depth.txt", emit: rarefaction_depth
    path "versions.yml", emit: versions

    script:
    """
    # Tạo thư mục kết quả
    mkdir -p iqtree_output

    # Kiểm tra phân bố số reads để chọn rarefaction depth (>80% mẫu)
    qiime feature-table summarize \\
        --i-table $feature_table \\
        --o-visualization feature-table.qzv
    
    # Tìm giá trị rarefaction depth tối ưu từ feature-table.qzv
    best_depth=\$(qiime tools peek feature-table.qzv | grep -A2 'Frequency per sample' | tail -n1 | awk '{print \$2 * 0.8}')
    best_depth=\$(printf "%.0f" "\$best_depth")  # Làm tròn số
    echo "Chọn rarefaction depth: \$best_depth" > rarefaction_depth.txt

    # Căn chỉnh lại trình tự đại diện trước khi tạo cây
    qiime phylogeny align-to-tree-mafft-fasttree \\
        --i-sequences $aligned_rep_seqs \\
        --o-alignment iqtree_output/aligned-sequences.qza \\
        --o-masked-alignment iqtree_output/masked-alignment.qza \\
        --o-tree iqtree_output/unrooted-tree.qza \\
        --o-rooted-tree iqtree_output/rooted-tree.qza

    mv iqtree_output/unrooted-tree.qza .
    mv iqtree_output/rooted-tree.qza .

    # Lưu thông tin phiên bản
    echo 'QIIME2_VERSION: ' \$(qiime --version) > versions.yml
    echo 'IQ-TREE_VERSION: ' \$(iqtree2 --version | head -n 1) >> versions.yml
    """
}
