#!/usr/bin/env Rscript
# Analysis of CTVT whole-genome variants
# Step 9.3. Process VEP annotation
# Adrian Baez-Ortega, 2018

# Load packages
library(stringr)

# Function to print numbers with commas
num = function(x) prettyNum(x, big.mark=",")


# Load variant data
load("VariantTables_Indices_Filt.RData")
snvs.metadata.host    <- readRDS("VariantTables_Split_snvs_Metadata.host.RDS")
indels.metadata.host  <- readRDS("VariantTables_Split_indels_Metadata.host.RDS")
snvs.metadata.tonly   <- readRDS("VariantTables_Split_Filt_snvs_Metadata.tonly.RDS")
indels.metadata.tonly <- readRDS("VariantTables_Split_Filt_indels_Metadata.tonly.RDS")
snvs.nv.tonly         <- readRDS("VariantTables_Split_Filt_snvs_NV.tonly.RDS")
indels.nv.tonly       <- readRDS("VariantTables_Split_Filt_indels_NV.tonly.RDS")

# Load annotation
annot.snvs.tonly = read.table("VEP_Annotation_SomaticSNVs_ToAnnotate.tsv", sep="\t",
                              header=T, skip=50, comment.char="", stringsAsFactors=F, check.names=F)
annot.snvs.host = read.table("VEP_Annotation_HostSNVs_ToAnnotate.tsv", sep="\t",
                             header=T, skip=50, comment.char="", stringsAsFactors=F, check.names=F)
annot.indels.tonly = read.table("VEP_Annotation_SomaticIndels_ToAnnotate.tsv", sep="\t",
                                header=T, skip=50, comment.char="", stringsAsFactors=F, check.names=F)
annot.indels.host = read.table("VEP_Annotation_HostIndels_ToAnnotate.tsv", sep="\t",
                               header=T, skip=50, comment.char="", stringsAsFactors=F, check.names=F)
stopifnot(identical(colnames(annot.snvs.tonly), colnames(annot.snvs.host)) &
              identical(colnames(annot.snvs.tonly), colnames(annot.indels.tonly)) &
              identical(colnames(annot.snvs.tonly), colnames(annot.indels.host)))

cat("Found annotation for:", num(length(unique(annot.snvs.tonly[,1]))),
    "/", num(nrow(snvs.metadata.tonly)), "tumour-only SNVs\n",
    "                    ", num(length(unique(annot.indels.tonly[,1]))),
    "/", num(nrow(indels.metadata.tonly)), "tumour-only indels\n",
    "                    ", num(length(unique(annot.snvs.host[,1]))),
    "/", num(nrow(snvs.metadata.host)), "host SNVs\n",
    "                    ", num(length(unique(annot.indels.host[,1]))),
    "/", num(nrow(indels.metadata.host)), "host indels\n")


# Process annotation
# Rename some of the columns
idx = match(c("#Uploaded_variation", "Allele", "Feature", "IMPACT", "DISTANCE", "STRAND", "SYMBOL",
              "SYMBOL_SOURCE", "BIOTYPE", "CANONICAL", "NEAREST", "EXON", "INTRON", "DOMAINS"),
            colnames(annot.snvs.tonly))
colnames(annot.snvs.tonly)[idx] =
    colnames(annot.snvs.host)[idx] =
    colnames(annot.indels.tonly)[idx] =
    colnames(annot.indels.host)[idx] = c("Meta_index", "Alt", "Transcript", "Impact", "Distance_to_tx",
                                         "Strand", "Symbol", "Symbol_source", "Biotype",
                                         "Canonical", "Nearest_TSS", "Exon", "Intron", "Domains")


# Discard uninformative fields
DISCARD = c("Feature_type", "CLIN_SIG", "SOMATIC", "PHENO", "MOTIF_NAME", "MOTIF_POS",
            "HIGH_INF_POS", "MOTIF_SCORE_CHANGE", "FLAGS", "GENE_PHENO")
cat("\nDiscarding uninformative fields:", paste(DISCARD, collapse=", "), "\n")

idx = -match(DISCARD, colnames(annot.snvs.tonly))
annot.snvs.tonly = annot.snvs.tonly[, idx]
annot.snvs.host = annot.snvs.host[, idx]
annot.indels.tonly = annot.indels.tonly[, idx]
annot.indels.host = annot.indels.host[, idx]


# For each gene: if annotated on multiple transcripts,
# select annotation on the canonical (or first) transcript
process.transcripts = function(annotation) {
    keep.idx = rep(FALSE, nrow(annotation))
    canonical.idx = annotation$Canonical == "YES"

    for (gene in unique(annotation$Gene)) {
        gene.idx = annotation$Gene == gene
        transcripts = unique(annotation$Transcript[gene.idx])

        # If there is only one transcript, keep it
        if (length(transcripts) == 1) {
            keep.idx[gene.idx] = TRUE
        }
        # If multiple transcripts, keep canonical or first annotated one
        else {
            gene.canon.idx = gene.idx & canonical.idx
            if (any(gene.canon.idx)) {
                if (length(unique(annotation$Transcript[gene.canon.idx])) > 1) {
                    stop("Multiple canonical transcripts for ", gene)
                }
                keep.idx[gene.canon.idx] = TRUE
            }
            else {
                cat("No canonical transcript for", gene, "- selecting first annotated transcript\n")
                keep.idx[annotation$Transcript == transcripts[1]] = TRUE
            }
        }
    }
    keep.idx
}

cat("\nSelecting canonical transcript annotation per variant:\n")
cat("Processing host SNVs...\n")
keep.snv.h.idx = process.transcripts(annot.snvs.host)
cat("Processing tumour-only SNVs...\n")
keep.snv.t.idx = process.transcripts(annot.snvs.tonly)
cat("Processing host indels...\n")
keep.ind.h.idx = process.transcripts(annot.indels.host)
cat("Processing tumour-only indels...\n")
keep.ind.t.idx = process.transcripts(annot.indels.tonly)

# Count lost variants
cat("\nVariants lost due to transcript selection:\n",
    num(sum(!(unique(annot.snvs.tonly$Meta_index) %in% unique(annot.snvs.tonly$Meta_index[keep.snv.t.idx])))),
    "tumour-only SNVs\n",
    num(sum(!(unique(annot.indels.tonly$Meta_index) %in% unique(annot.indels.tonly$Meta_index[keep.ind.t.idx])))),
    "tumour-only indels\n",
    num(sum(!(unique(annot.snvs.host$Meta_index) %in% unique(annot.snvs.host$Meta_index[keep.snv.h.idx])))),
    "host SNVs\n",
    num(sum(!(unique(annot.indels.host$Meta_index) %in% unique(annot.indels.host$Meta_index[keep.ind.h.idx])))),
    "host indels\n")

# Discard unselected transcripts
annot.snvs.tonly = annot.snvs.tonly[keep.snv.t.idx, ]
annot.snvs.host = annot.snvs.host[keep.snv.h.idx, ]
annot.indels.tonly = annot.indels.tonly[keep.ind.t.idx, ]
annot.indels.host = annot.indels.host[keep.ind.h.idx, ]


# Add extra fields
cat("\nAdding extra fields:\n")
cat("Chrom\n")
annot.snvs.tonly$Chrom = snvs.metadata.tonly$CHROM[annot.snvs.tonly$Meta_index]
annot.snvs.host$Chrom = snvs.metadata.host$CHROM[annot.snvs.host$Meta_index]
annot.indels.tonly$Chrom = indels.metadata.tonly$CHROM[annot.indels.tonly$Meta_index]
annot.indels.host$Chrom = indels.metadata.host$CHROM[annot.indels.host$Meta_index]

cat("Pos\n")  # adjusted for indels
annot.snvs.tonly$Position = snvs.metadata.tonly$POS[annot.snvs.tonly$Meta_index]
annot.snvs.host$Position = snvs.metadata.host$POS[annot.snvs.host$Meta_index]
annot.indels.tonly$Position = indels.metadata.tonly$POS[annot.indels.tonly$Meta_index] + 1
annot.indels.host$Position = indels.metadata.host$POS[annot.indels.host$Meta_index] + 1

cat("Ref\n")  # adjusted for indels
annot.snvs.tonly$Ref = snvs.metadata.tonly$REF[annot.snvs.tonly$Meta_index]
annot.snvs.host$Ref = snvs.metadata.host$REF[annot.snvs.host$Meta_index]
annot.indels.tonly$Ref = substring(indels.metadata.tonly$REF[annot.indels.tonly$Meta_index], 2)
annot.indels.host$Ref = substring(indels.metadata.host$REF[annot.indels.host$Meta_index], 2)
annot.indels.tonly$Ref[annot.indels.tonly$Ref == ""] = "-"
annot.indels.host$Ref[annot.indels.host$Ref == ""] = "-"

# cat("CDS_Ref, CDS_Alt\n")  # only SNVs
# annot.indels.tonly$CDS_Ref = annot.indels.tonly$CDS_Alt =
#     annot.indels.host$CDS_Ref = annot.indels.host$CDS_Alt =
#     annot.snvs.tonly$CDS_Ref = annot.snvs.tonly$CDS_Alt =
#     annot.snvs.host$CDS_Ref = annot.snvs.host$CDS_Alt = "-"
# for (i in 1:nrow(annot.snvs.tonly)) {
#     codon = annot.snvs.tonly$Codons[i]
#     if (codon != "-") {
#         bases = str_split_fixed(codon, "", 7)
#         bases = bases[grep("[[:upper:]]", bases)]
#         annot.snvs.tonly$CDS_Ref[i] = bases[1]
#         annot.snvs.tonly$CDS_Alt[i] = bases[2]
#     }
# }
# for (i in 1:nrow(annot.snvs.host)) {
#     codon = annot.snvs.host$Codons[i]
#     if (codon != "-") {
#         bases = str_split_fixed(codon, "", 7)
#         bases = bases[grep("[[:upper:]]", bases)]
#         annot.snvs.host$CDS_Ref[i] = bases[1]
#         annot.snvs.host$CDS_Alt[i] = bases[2]
#     }
# }

cat("Ref_trinucleotide\n")  # only SNVs
annot.indels.tonly$Ref_trinucleotide = "-"
annot.indels.host$Ref_trinucleotide = "-"
annot.snvs.tonly$Ref_trinucleotide = substr(str_split_fixed(snvs.metadata.tonly$INFO[annot.snvs.tonly$Meta_index],
                                                            ";", 20)[, 12],
                                            13, 15)
annot.snvs.host$Ref_trinucleotide = substr(str_split_fixed(snvs.metadata.host$INFO[annot.snvs.host$Meta_index],
                                                            ";", 20)[, 12],
                                            13, 15)

cat("Length\n")  # only indels; sign indicates ins/del
annot.snvs.tonly$Length = "-"
annot.snvs.host$Length = "-"
annot.indels.tonly$Length = nchar(indels.metadata.tonly$ALT[annot.indels.tonly$Meta_index]) -
    nchar(indels.metadata.tonly$REF[annot.indels.tonly$Meta_index])
annot.indels.host$Length = nchar(indels.metadata.host$ALT[annot.indels.host$Meta_index]) -
    nchar(indels.metadata.host$REF[annot.indels.host$Meta_index])

cat("Number_tumours\n")
annot.snvs.host$Number_tumours = "-"
annot.indels.host$Number_tumours = "-"
annot.snvs.tonly$Number_tumours = rowSums(snvs.nv.tonly[annot.snvs.tonly$Meta_index, tumours] >= MIN.READS)
annot.indels.tonly$Number_tumours = rowSums(indels.nv.tonly[annot.indels.tonly$Meta_index, tumours] >= MIN.READS)


# Rearrange columns
stopifnot(identical(colnames(annot.snvs.tonly), colnames(annot.snvs.host)) &
              identical(colnames(annot.snvs.tonly), colnames(annot.indels.tonly)) &
              identical(colnames(annot.snvs.tonly), colnames(annot.indels.host)))
idx = match(c("Meta_index", "Location", "Chrom", "Position", "Ref", "Alt", "Length",  # "CDS_Ref", "CDS_Alt",
              "Ref_trinucleotide", "Number_tumours", "Gene", "Transcript", "Symbol", "cDNA_position",
              "CDS_position", "Protein_position", "Exon", "Intron", "Strand", "Amino_acids", "Codons",
              "Consequence", "Impact", "SIFT", "Domains", "Nearest_TSS", "Distance_to_tx", "Biotype", "Canonical",
              "Existing_variation", "HGNC_ID", "Symbol_source"),
            colnames(annot.snvs.tonly))

annot.snvs.tonly = annot.snvs.tonly[, idx]
annot.snvs.host = annot.snvs.host[, idx]
annot.indels.tonly = annot.indels.tonly[, idx]
annot.indels.host = annot.indels.host[, idx]


# Save tables
cat("\nSaving processed annotation tables to data/ directory...\n")
dir.create("data")
saveRDS(annot.snvs.tonly,   file = "data/AnnotationTables_annot.snvs.tonly.RDS")
saveRDS(annot.snvs.host,    file = "data/AnnotationTables_annot.snvs.host.RDS")
saveRDS(annot.indels.tonly, file = "data/AnnotationTables_annot.indels.tonly.RDS")
saveRDS(annot.indels.host,  file = "data/AnnotationTables_annot.indels.host.RDS")
