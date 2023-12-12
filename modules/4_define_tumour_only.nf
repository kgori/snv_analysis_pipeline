process define_tumour_only {
    input:
    path rds
    path samplenames
    path matchedhosts
    path purities

    output:
    path "data/*.RData"
    path "data/*.RDS"

    script:
    """
    4_define_tumour_only.R ${samplenames} ${matchedhosts} ${purities}
    """
}
