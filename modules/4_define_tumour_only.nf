process define_tumour_only {
    input:
    path rds
    path samplenames
    path matchedhosts
    path purities

    output:
    path "data/*.RDS"
    path "data/*.RData"

    script:
    """
    define_tumour_only.R ${samplenames} ${matchedhosts} ${purities}
    """
}
