params.SNVsFile = ""
params.indelsFile = ""
params.sampleList = ""
params.panelFile = ""
params.outputDir = "./out"

include { vcf_to_read_counts }     from "./modules/snv_postprocessing.nf" 
include { assign_germline_status } from "./modules/snv_postprocessing.nf"
include { extract_somatic }        from "./modules/build_guide_tree.nf"
include { make_fasta }             from "./modules/build_guide_tree.nf"
include { build_guide_tree }       from "./modules/build_guide_tree.nf"
include { extract_vcf_data }       from "./modules/1_extract_info.nf"
include { import_variants }        from "./modules/2_import_variants.nf"
include { filter_variants }        from "./modules/3_filter_variants.nf"

workflow {
    println "Variants File = ${params.SNVsFile}"
    println "Indels File = ${params.indelsFile}"
    println "Sample List = ${params.sampleList}"
    println "Panel File = ${params.panelFile}"
    println "Output Dir = ${params.outputDir}"
    variants = Channel.fromPath("${params.SNVsFile}", checkIfExists: true)
    indels = Channel.fromPath("${params.indelsFile}", checkIfExists: true)
    samples = Channel.fromPath("${params.sampleList}", checkIfExists: true)
    panel = Channel.fromPath("${params.panelFile}", checkIfExists: true)

    counts = vcf_to_read_counts(variants)
    annotated_variants = assign_germline_status(counts, samples, panel)

    somatic_snvs = extract_somatic(annotated_variants[0])
    somatic_fasta = make_fasta(somatic_snvs, samples)
    guide_tree = build_guide_tree(somatic_fasta)

    extract_ch =
        variants.map { it -> tuple("snvs", it) }
        .concat(indels.map { it -> tuple("indels", it) })
    extract = extract_vcf_data(extract_ch)
    imported = extract | flatten | import_variants
    filtered = filter_variants(imported | collect, samples, guide_tree)
}
