process BUILD_PHYLOGENETIC_TREE {
    tag "qiime2_phylogenetic_tree"
    container "quay.io/qiime2/amplicon:2024.5"

    input:
    path rep_seqs  // 🔹 Nhận đầu vào từ FILTER_SEQS

    output:
    path "aligned_rep_seqs.qza", emit: aligned_rep_seqs
    path "masked_aligned_rep_seqs.qza", emit: masked_aligned_rep_seqs
    path "unrooted_tree.qza", emit: unrooted_tree
    path "rooted_tree.qza", emit: rooted_tree
    publishDir "qiime_out", mode: 'copy'

    script:
    """
    # 🔹 Căn chỉnh trình tự bằng MAFFT
    qiime alignment mafft \\
        --i-sequences ${rep_seqs} \\
        --o-alignment aligned_rep_seqs.qza

    # 🔹 Lọc nhiễu bằng Mask
    qiime alignment mask \\
        --i-alignment aligned_rep_seqs.qza \\
        --o-masked-alignment masked_aligned_rep_seqs.qza

    # 🔹 Xây dựng cây phát sinh loài chưa có gốc
    qiime phylogeny fasttree \\
        --i-alignment masked_aligned_rep_seqs.qza \\
        --o-tree unrooted_tree.qza

    # 🔹 Định gốc cây
    qiime phylogeny midpoint-root \\
        --i-tree unrooted_tree.qza \\
        --o-rooted-tree rooted_tree.qza
    qiime tools export \
        --input-path rooted_tree.qza \
        --output-path exported_tree

    """
}

