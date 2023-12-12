process filter_tumour_only {
    input:
    path indices
    path rds
    path host_filter_list
    path repeats
    path purity
    path matchedhosts

    output:
    path "data/VariantTables_Indices_Filt.RData"
    path "data/*.RDS"

    script:
    """
    6_FilterTOnly.R ${matchedhosts}
    """
}
