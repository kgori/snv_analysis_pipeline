process filter_variants {
    input:
    path rds
    path samplenames
    path tree

    output:
    path("data/*.RDS")
    path("site_patterns/*.pdf")

    publishDir "${params.outputDir}/site_patterns", pattern: 'site_patterns/*.pdf', mode: 'copy'

    script:
    """
    filter_variants.R ${samplenames} ${tree}
    """
}
