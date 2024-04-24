process prep_vcf {
    input:
    path rds

    output:
    path "data/*.vcf"

    script:
    """
    9.1_AnnotateVariants_PrepVCF.R
    """
}

/* NB: You need to download an ensembl cache from 
   https://ftp.ensembl.org/pub/release-104/variation/indexed_vep_cache/
   (or whatever the latest release is), and set its location in the
   nextflow.config's containerOptions string as '-B path-to-cache:/cache
   TODO: make this a parameter, not something hard-coded in a config file
*/
process run_vep {
    input:
    path vcf

    output:
    path "data/*.tsv"

    publishDir "${params.outputDir}/vep_output/raw", mode: 'copy'

    script:
    """
    mkdir -p data
    9.2_RunVEP.sh ${task.cpus} data ${vcf}
    """
}

process process_vep_annotation {
    input:
    path filtered_index
    path filtered_rds
    path tonly_rds
    path vep_tsvs

    output:
    path "data/*.RDS"

    publishDir "${params.outputDir}/vep_output/processed", mode: 'copy'

    script:
    """
    9.3_ProcessAnnotation.R
    """
} 
