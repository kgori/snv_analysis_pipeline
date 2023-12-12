params.SNVsFile = ""
params.indelsFile = ""
params.sampleList = ""
params.hostMatch = ""
params.panelFile = ""
params.outputDir = "./out"
params.reference = "" /* Reference genome */
params.lowcovDir = "" /* Folder containing bams or crams of low coverage sequencing */
params.repeatsFile = "" /* RData file of simple repeat regions for the reference */
params.trinucsFile = "" /* RData file of trinucleotide contexts for the reference */

def remove_duplicate_filepair_keys(primary_ch, secondary_ch) {
    // primary_ch and secondary_ch are channels produced by
    // `fromFilePairs`. If any keys are present in both channels,
    // remove the key from the secondary channel.
    // This will avoid duplicating work on a bam file if it has
    // a bai and a csi index.
    return primary_ch.concat(secondary_ch).groupTuple() |
        map { it -> tuple(it[0], it[1][0]) }
}


include { vcf_to_read_counts }      from "./modules/0_snv_postprocessing.nf" 
include { assign_germline_status }  from "./modules/0_snv_postprocessing.nf"
include { extract_somatic }         from "./modules/0_build_guide_tree.nf"
include { make_fasta }              from "./modules/0_build_guide_tree.nf"
include { build_guide_tree }        from "./modules/0_build_guide_tree.nf"
include { extract_vcf_data }        from "./modules/1_extract_info.nf"
include { import_variants }         from "./modules/2_import_variants.nf"
include { filter_variants }         from "./modules/3_filter_variants.nf"
include { define_tumour_only }      from "./modules/4_define_tumour_only.nf"
include { prepare_vcf }             from "./modules/5_genotype_host_panel.nf"
include { genotype_low_cov_hosts }  from "./modules/5_genotype_host_panel.nf"
include { make_filter_list }        from "./modules/5_genotype_host_panel.nf"
include { filter_tumour_only }      from "./modules/6_filter_tumour_only"
include { vafhist_spectra }         from "./modules/7_vafhist_spectra.nf"
/* Skipping section 8 as it is only for troubleshooting */
include { prep_vcf }                from "./modules/9_vep.nf"
include { run_vep }                 from "./modules/9_vep.nf"
include { process_vep_annotation }  from "./modules/9_vep.nf"
include { make_alignment }          from "./modules/10_phylotree.nf"
include { run_unpartitioned_raxml } from "./modules/10_phylotree.nf"
include { run_partitioned_raxml }   from "./modules/10_phylotree.nf"

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
    println "Trinucleotide Contexts File = ${params.trinucsFile}"

    variants = Channel.fromPath("${params.SNVsFile}", checkIfExists: true)
    indels = Channel.fromPath("${params.indelsFile}", checkIfExists: true)
    samples = Channel.fromPath("${params.sampleList}", checkIfExists: true)
    hostmatch = Channel.fromPath("${params.hostMatch}", checkIfExists: true)
    panel = Channel.fromPath("${params.panelFile}", checkIfExists: true)
    repeats = Channel.fromPath("${params.repeatsFile}", checkIfExists: true)
    trinucs = Channel.fromPath("${params.trinucsFile}", checkIfExists: true)
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
    tonly_index = tonly[0]
    tonly_rds = tonly[1]
    
    prepvcf = prepare_vcf(tonly_rds)

    // Genotyping
    // geno = genotype_low_cov_hosts(low_cov_files.combine(reference).combine(prepvcf))
    // filter_list = make_filter_list(geno.collect())

    // Filtering
    filtered_tumours = filter_tumour_only(tonly_index,
        tonly_rds.collect(), panel, repeats, purities, hostmatch) 
    filt_index = filtered_tumours[0]
    filt_rds = filtered_tumours[1]

    // VAF Histograms
    spectra = vafhist_spectra(filt_index,
        tonly_rds.collect(), 
        filt_rds.collect(),
        trinucs,
        hostmatch)

    spectra_tsvs = spectra[0]
    spectra_plots = spectra[1]
    spectra_data = spectra[2]
    spectra_purities = spectra[3]

    // VEP
    vep_prep = prep_vcf(tonly_rds.collect().combine(filt_rds.collect()))
    vep_result = run_vep(vep_prep)
    vep_processed = process_vep_annotation(filt_index, filt_rds.collect(), tonly_rds.collect(), vep_result)

    // Phylogeny
    alignment = make_alignment(spectra_purities, filt_index, filt_rds.collect(), vep_processed)
    phy = alignment[0]
    partition = alignment[2]
    run_unpartitioned_raxml(phy)
    run_partitioned_raxml(phy, partition)
}
