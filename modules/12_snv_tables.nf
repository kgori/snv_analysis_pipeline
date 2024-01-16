process annotate_metadata {
    input:
    path metadata
    path host_metadata
    path tumour_metadata

    output:
    path "data/VariantTables_Filtered_snvs_Annotated_Metadata.RDS", emit: snvs
    path "data/VariantTables_Filtered_indels_Annotated_Metadata.RDS", emit: indels

    script:
    """
    12.1_Annotate_Context_Metadata.R
    """
}

process write_arrow {
    input:
    path annotated_snv_metadata
    path annotated_indel_metadata
    path data
    path matchedhosts

    output:
    path "data/dataset"

    script:
    """
    12.2_Write_Arrow.R ${matchedhosts}
    """
}
