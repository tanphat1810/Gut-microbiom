process BETA_DIVERSITY {
    tag "qiime2_beta_diversity"

    conda "bioconda::qiime2=2023.7"
    container "quay.io/qiime2/amplicon:2024.5"

    input:
    path feature_table  // Feature table đầu vào
    path rooted_tree    // Cây phát sinh loài (cần cho UniFrac)

    output:
    path "jaccard_distance_matrix.qza", emit: jaccard
    path "braycurtis_distance_matrix.qza", emit: bray_curtis
    path "unweighted_unifrac_distance_matrix.qza", emit: unweighted_unifrac
    path "weighted_unifrac_distance_matrix.qza", emit: weighted_unifrac
    publishDir "qiime_out", mode: 'copy'

    script:
    """
    # 🔹 Tính toán Jaccard & Bray-Curtis
    for metric in 'jaccard' 'braycurtis'; do
        qiime diversity beta \\
            --i-table ${feature_table} \\
            --p-metric \$metric \\
            --o-distance-matrix \${metric}_distance_matrix.qza
    done

    # 🔹 Tính toán UniFrac
    for metric in 'unweighted_unifrac' 'weighted_unifrac'; do
        qiime diversity beta-phylogenetic \\
            --i-table ${feature_table} \\
            --i-phylogeny ${rooted_tree} \\
            --p-metric \$metric \\
            --o-distance-matrix \${metric}_distance_matrix.qza
    done
    """
}

