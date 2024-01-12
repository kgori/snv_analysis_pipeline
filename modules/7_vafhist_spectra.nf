process vafhist_spectra {
    input:
    path indices
    path rds
    path filtered_rds
    path trinucs
    path matchedhosts

    output:
    path "output/*.tsv"
    path "output/*.pdf"
    path "data/VariantSpectra_Filt.RData"
    path "data/TumourPurity.RData"

    publishDir "${params.outputDir}/VAFhist_Spectra", mode: 'copy'

    script:
    """
    7_VAFhist_Spectra.R ${matchedhosts}
    """
}
