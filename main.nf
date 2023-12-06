params.SNVsFile = ""
params.indelsFile = ""
params.sampleList = ""
params.hostMatch = ""
params.panelFile = ""
params.outputDir = "./out"
params.reference = "" /* Reference genome */
params.lowcovDir = "" /* Folder containing bams or crams of low coverage sequencing */
params.repeatsFile = "" /* RData file of simple repeat regions for the reference */

def remove_duplicate_filepair_keys(primary_ch, secondary_ch) {
    // primary_ch and secondary_ch are channels produced by
    // `fromFilePairs`. If any keys are present in both channels,
    // remove the key from the secondary channel.
    // This will avoid duplicating work on a bam file if it has
    // a bai and a csi index.
    return primary_ch.concat(secondary_ch).groupTuple() |
        map { it -> tuple(it[0], it[1][0]) }
}


include { vcf_to_read_counts }     from "./modules/snv_postprocessing.nf" 
include { assign_germline_status } from "./modules/snv_postprocessing.nf"
include { extract_somatic }        from "./modules/build_guide_tree.nf"
include { make_fasta }             from "./modules/build_guide_tree.nf"
include { build_guide_tree }       from "./modules/build_guide_tree.nf"
include { extract_vcf_data }       from "./modules/1_extract_info.nf"
include { import_variants }        from "./modules/2_import_variants.nf"
include { filter_variants }        from "./modules/3_filter_variants.nf"
include { define_tumour_only }     from "./modules/4_define_tumour_only.nf"
include { prepare_vcf }            from "./modules/5_genotype_host_panel.nf"
include { genotype_low_cov_hosts } from "./modules/5_genotype_host_panel.nf"
include { make_filter_list }       from "./modules/5_genotype_host_panel.nf"
include { filter_tumour_only }     from "./modules/6_filter_tumour_only"

workflow {
    println "Variants File = ${params.SNVsFile}"
    println "Indels File = ${params.indelsFile}"
    println "Sample List = ${params.sampleList}"
    println "Matched Hosts = ${params.hostMatch}"
    println "Panel File = ${params.panelFile}"
    println "Output Dir = ${params.outputDir}"
    println "Reference = ${params.reference}"
    println "Low Coverage Dir = ${params.lowCovDir}"
    println "Repeats File = ${params.repeatsFile}"

    variants = Channel.fromPath("${params.SNVsFile}", checkIfExists: true)
    indels = Channel.fromPath("${params.indelsFile}", checkIfExists: true)
    samples = Channel.fromPath("${params.sampleList}", checkIfExists: true)
    hostmatch = Channel.fromPath("${params.hostMatch}", checkIfExists: true)
    panel = Channel.fromPath("${params.panelFile}", checkIfExists: true)
    repeats = Channel.fromPath("${params.repeatsFile}", checkIfExists: true)
    bam_bai_inputs   = Channel.fromFilePairs("${params.lowCovDir}/*.bam{,.bai}")
    bam_csi_inputs   = Channel.fromFilePairs("${params.lowCovDir}/*.bam{,.csi}")
    bam_inputs = remove_duplicate_filepair_keys(bam_csi_inputs, bam_bai_inputs)
    cram_inputs = Channel.fromFilePairs("${params.lowCovDir}/*.cram{,.crai}")
    low_cov_files = bam_inputs.concat(cram_inputs)
    
    refGlob = "${params.reference}" + "{,.fai}"
    reference = Channel.fromFilePairs("${refGlob}", checkIfExists: true)

    counts = vcf_to_read_counts(variants)
    annotated_variants = assign_germline_status(counts, samples, panel)
    purities = annotated_variants[2]

    somatic_snvs = extract_somatic(annotated_variants[0])
    somatic_fasta = make_fasta(somatic_snvs, samples)
    guide_tree = build_guide_tree(somatic_fasta)

    extract_ch =
        variants.map { it -> tuple("snvs", it) }
        .concat(indels.map { it -> tuple("indels", it) })
    extract = extract_vcf_data(extract_ch)
    imported = extract | flatten | import_variants
    filtered = filter_variants(imported | collect, samples, guide_tree)

    tonly = define_tumour_only(filtered[0], samples, hostmatch, purities)

    prepvcf = prepare_vcf(tonly[0])

    // Genotyping
    geno = genotype_low_cov_hosts(low_cov_files.combine(reference).combine(prepvcf))
    filter_list = make_filter_list(geno.collect())

    // Filtering
    filtered_tumours = filter_tumour_only(tonly[0],
        tonly[1].collect(), filter_list, repeats, purities, hostmatch) 

}
