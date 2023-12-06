process prepare_vcf {
    input:
    path vcf

    output:
    tuple val("genotyping_vcf"), path("data/SomaticVariants_ToGenotype.vcf.gz*")

    script:
    """
    5.1_GenotypeHostPanel_PrepVCF.R
    
    sort -k1,1 -k2,2n -k4,4 -k5,5 data/SomaticVariants_ToGenotype.vcf \
      > data/SomaticVariants_ToGenotype.sorted.vcf \
      && mv data/SomaticVariants_ToGenotype.sorted.vcf data/SomaticVariants_ToGenotype.vcf
    
    cat <(echo '##fileformat=VCFv4.1\n##source=Somatypus\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO')\
      data/SomaticVariants_ToGenotype.vcf\
      | bgzip > data/SomaticVariants_ToGenotype.vcf.gz \
      && rm data/SomaticVariants_ToGenotype.vcf
   
    tabix -C data/SomaticVariants_ToGenotype.vcf.gz
    """
}

process genotype_low_cov_hosts {
    input:
    tuple val(bamID), path(alignmentFiles),
      val(refID), path (referenceFiles),
      val(vcfID), path(source_vcf)

    output:
    path "genotyped/${bamID}_FINAL.vcf.gz"

    script:
    """
    5.2_GenotypeVariantsVCF.sh ${source_vcf[0]} . genotyped ${bamID} ${task.cpus} ${referenceFiles[0]}
    """
}

process make_filter_list {
    input:
    path vcf

    output:
    path "filter_list.txt"

    script:
    """
    make_filter_list.py ${vcf} > log.txt
    """
}
