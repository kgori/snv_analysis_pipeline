process extract_vcf_data {
    input:
    path vcf

    output:
    tuple path("${vcf.getSimpleName()}_Metadata.tsv.gz"), \
      path("${vcf.getSimpleName()}_NR.tsv.gz"), \
      path("${vcf.getSimpleName()}_NV.tsv.gz")

    script:
    """
    extract_vcf_data.py ${vcf} .
    """
}
