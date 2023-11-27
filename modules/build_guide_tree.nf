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
    path samplelist

    output:
    path "somatic_variants.fa"

    script:
    """
    somatic_fasta_from_tsv.py ${snvs} ${samplelist} somatic_variants.fa
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
    iqtree -T ${task.cpus} -s ${fasta}
    """
}
