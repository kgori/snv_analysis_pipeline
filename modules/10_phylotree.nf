process make_alignment {
    input:
    path purity
    path filtered_index
    path filtered_rds
    path annotated_rds

    output:
    path "data/TumourAlignment.phylip"
    path "data/TumourAlignment_LastCodingSite.txt"
    path "data/partitions.txt"

    script:
    """
    10.1_BuildAlignment.R
    """
}

process run_unpartitioned_raxml {
    label 'raxml'

    input:
    path phy

    output:
    path "data/*"

    publishDir "${params.outputDir}/raxml_unpartitioned"

    script:
    """
    10.2_RunRAxML.sh ${phy} \$(realpath data) ${task.cpus} no_partition
    """
}

process run_partitioned_raxml {
    label 'raxml'

    input:
    path phy
    path partitions

    output:
    path "data/*"

    publishDir "${params.outputDir}/raxml_partitioned"

    script:
    """
    10.2_RunRAxML.sh ${phy} \$(realpath data) ${task.cpus} ${partitions}
    """
}
