params.variantsFile = ""
params.sampleList = ""
params.panelFile = ""
params.outputDir = "./out"

include { vcf_to_read_counts }     from "./modules/snv_postprocessing.nf" 
include { assign_germline_status } from "./modules/snv_postprocessing.nf"
include { extract_somatic }        from "./modules/build_guide_tree.nf"
include { make_fasta }             from "./modules/build_guide_tree.nf"

workflow {
    println "Variants File = ${params.variantsFile}"
    println "Sample List = ${params.sampleList}"
    println "Panel File = ${params.panelFile}"
    println "Output Dir = ${params.outputDir}"
    variants = Channel.fromPath("${params.variantsFile}", checkIfExists: true)
    samples = Channel.fromPath("${params.sampleList}", checkIfExists: true)
    panel = Channel.fromPath("${params.panelFile}", checkIfExists: true)

    counts = vcf_to_read_counts(variants)
    annotated_variants = assign_germline_status(counts, samples, panel)

    somatic_snvs = extract_somatic(annotated_variants[0])
    somatic_fasta = make_fasta(somatic_snvs, samples)
}
