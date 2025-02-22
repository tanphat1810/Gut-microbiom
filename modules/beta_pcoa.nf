process BETA_PCOA {
    tag "qiime2_beta_pcoa"

    conda "bioconda::qiime2=2023.7"
    container "quay.io/qiime2/amplicon:2024.5"

    input:
    path jaccard
    path bray_curtis
    path unweighted_unifrac
    path weighted_unifrac

    output:
    path "*_pcoa.qza", emit: pcoa_results
    publishDir "qiime_out", mode: 'copy'

 script:
"""
qiime diversity pcoa \\
    --i-distance-matrix ${jaccard} \\
    --o-pcoa jaccard_pcoa.qza

qiime diversity pcoa \\
    --i-distance-matrix ${bray_curtis} \\
    --o-pcoa bray_curtis_pcoa.qza

qiime diversity pcoa \\
    --i-distance-matrix ${unweighted_unifrac} \\
    --o-pcoa unweighted_unifrac_pcoa.qza

qiime diversity pcoa \\
    --i-distance-matrix ${weighted_unifrac} \\
    --o-pcoa weighted_unifrac_pcoa.qza
"""

   
}

