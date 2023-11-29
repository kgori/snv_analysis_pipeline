process filter_variants {
    input:
    path rds
    path samplenames
    path tree

    output:
    path("*filtered.RDS")

    script:
    """
    filter_variants.R ${samplenames} ${tree}
    """
}
