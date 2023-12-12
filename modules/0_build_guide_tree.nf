process extract_somatic {
    input:
    path snvs

    output:
    path "somatic_snvs.tsv.gz"

    script:
    """
    gunzip -c ${snvs} | grep "CHROM\\|somatic" | rev | cut -d'\t' -f2- | rev | gzip -c > somatic_snvs.tsv.gz
    """
}

process make_fasta {
    input:
    path snvs
    path samplenames

    output:
    path "somatic_variants.fa"

    script:
    """
    0.2_somatic_fasta_from_tsv.py ${snvs} ${samplenames} somatic_variants.fa
    """
}

process build_guide_tree {
    input:
    path fasta

    output:
    path "somatic_variants.fa.treefile"

    publishDir "${params.outputDir}/guide_tree", mode: 'copy'

    script:
    """
    iqtree2 -T ${task.cpus} -m GTR+G{4} -s ${fasta} 
    """
}
