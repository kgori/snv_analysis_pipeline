process import_variants {
    input:
    path tsv

    output:
    path "VariantTables_${tsv.getSimpleName()}.RData"

    script:
    """
    import_variants.R ${tsv} VariantTables_${tsv.getSimpleName()}.RData
    """
}
