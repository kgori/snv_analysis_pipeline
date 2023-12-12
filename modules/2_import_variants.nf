process import_variants {
    input:
    path tsv

    output:
    path "VariantTables_${tsv.getSimpleName()}.RDS"

    script:
    """
    2_import_variants.R ${tsv} VariantTables_${tsv.getSimpleName()}.RDS
    """
}
