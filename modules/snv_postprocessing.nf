process vcf_to_read_counts {
    input:
    path vcf

    output:
    path "snvs.tsv.gz"

    shell:
    """
    bcftools query -l ${vcf} 2>/dev/null | \
        tr '\\n' '\\t' | \
        awk 'BEGIN{printf("#CHROM\\tPOS\\tREF\\tALT\\tQUAL\\tINFO");}\
             {for(i=1; i<=NF; i++) printf("\\t%s.nr", \$i);}\
             {for (i=1; i<=NF; i++) printf("\\t%s.nv", \$i);}\
             END{printf("\\n");}' | gzip -c > snvs.tsv.gz
    bcftools query -f '%CHROM\\t%POS\\t%REF\\t%ALT\\t%QUAL\\t%SC[\\t%NR][\\t%NV]\\n' ${vcf} | gzip -c >> snvs.tsv.gz 
    """
}

process assign_germline_status {
    input:
    path tsv
    path samplenames
    path panel

    output:
    path "snvs_status.tsv.gz"
    path "snvs_status.tsv.gz.pdf"
    path "TumourPurity.RData"

    publishDir "${params.outputDir}/snvs_status"

    shell:
    """
    assign_germline.R ${tsv} ${samplenames} snvs_status.tsv.gz ${panel}
    """
}
