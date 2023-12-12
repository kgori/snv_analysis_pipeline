process extract_vcf_data {
    input:
    tuple val(id), path(vcf)

    output:
    tuple path("${id}_Metadata.tsv.gz"), \
      path("${id}_NR.tsv.gz"), \
      path("${id}_NV.tsv.gz")

    script:
    """
    1_extract_vcf_data.py ${vcf} ${id} .
    """
}
