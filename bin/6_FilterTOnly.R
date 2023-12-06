#!/usr/bin/env Rscript
# Analysis of CTVT whole-genome variants

# Step 6. Additional filtering of tumour-only variants
# Adrian Baez-Ortega, 2018

# Filters:
# 1. Low support filter: Discard variants with <5 supporting reads in all tumours.
# 2. Host support filter: Discard variants with >4 supporting reads in ≥1 host.
# 3. Mismapping filter: Discard variants with 1-2 supporting reads in >10 tumours.
# 4. Host panel filter: Discard variants with >4 total supporting reads across 
#    the panel of low-coverage hosts.
# 5. Repeat filter: Discard variants in, or within 5 bp of, annotated repeat regions
#    (we consider simple, low-complexity and tandem repeats, excluding exons).
# 6. Strand bias filter: Discard variants with ≤20% of total supporting reads on 
#    either strand.
# 7. Unmatched tumour-unique filter: Discard variants which are unique to an
#    unmatched tumour, and have VAF < 0.3 * purity.
#    NOTE: This filter can only be applied after obtaining purity estimates in
#    Step 7 (7_VAFhist_Spectra.R), and should be commented out when running this
#    script for the first time.


# Load packages
library(stringr)
library(data.table)
# library(sigfit, quietly=T)
suppressPackageStartupMessages(library(GenomicRanges))

# Function to print numbers with commas
num = function(x) prettyNum(x, big.mark=",")

# Load somatic variant data
load("VariantTables_Indices.RData")
snvs.metadata.tonly   <- readRDS("VariantTables_Split_snvs_Metadata.tonly.RDS")
snvs.nr.tonly         <- readRDS("VariantTables_Split_snvs_NR.tonly.RDS")
snvs.nv.tonly         <- readRDS("VariantTables_Split_snvs_NV.tonly.RDS")
snvs.vaf.tonly        <- readRDS("VariantTables_Split_snvs_VAF.tonly.RDS")
indels.metadata.tonly <- readRDS("VariantTables_Split_indels_Metadata.tonly.RDS")
indels.nr.tonly       <- readRDS("VariantTables_Split_indels_NR.tonly.RDS")
indels.nv.tonly       <- readRDS("VariantTables_Split_indels_NV.tonly.RDS")
indels.vaf.tonly      <- readRDS("VariantTables_Split_indels_VAF.tonly.RDS")

# Load host panel read counts
panel.nv.total = read.table("filter_list.txt", 
                            header=T, stringsAsFactors=F)
cat(num(nrow(panel.nv.total)), "variants genotyped on host panel\n\n")


# Prepare filter indices
cat("Before filtering:", num(nrow(snvs.metadata.tonly)), "tumour-only SNVs and",
    num(nrow(indels.metadata.tonly)), "tumour-only indels\n\n")

# 1) Low support filter
# Discard variants with <5 supporting reads in all tumours
# (i.e. keep variants with ≥5 reads in ≥1 tumour)
MIN.NV = 5
five.reads.snv.idx = apply(snvs.nv.tonly[, tumours], 1, max) >= MIN.NV
five.reads.ind.idx = apply(indels.nv.tonly[, tumours], 1, max) >= MIN.NV


# 2) Host support filter
# Discard variants with >4 supporting reads in ≥1 host
# (i.e. keep variants with ≤4 reads in all hosts)
MAX.NV = 4
no.host.supp.snv.idx = apply(snvs.nv.tonly[, hosts], 1, max) <= MAX.NV
no.host.supp.ind.idx = apply(indels.nv.tonly[, hosts], 1, max) <= MAX.NV


# 3) Mismapping filter
# Discard variants with 1-2 supporting reads in >10 tumours
MAX.NV = 2
MAX.TUMOURS = 10
not.mismap.snv.idx = rowSums(snvs.nv.tonly[, tumours] > 0 & 
                                    snvs.nv.tonly[, tumours] <= MAX.NV) <= MAX.TUMOURS
not.mismap.ind.idx = rowSums(indels.nv.tonly[, tumours] > 0 & 
                                    indels.nv.tonly[, tumours] <= MAX.NV) <= MAX.TUMOURS


# 4) Host panel filter
# Discard variants with >4 total supporting reads across the panel of low-coverage hosts
snv.str = paste(snvs.metadata.tonly$CHROM, snvs.metadata.tonly$POS,
                snvs.metadata.tonly$REF, snvs.metadata.tonly$ALT,
                sep=":")
ind.str = paste(indels.metadata.tonly$CHROM, indels.metadata.tonly$POS,
                indels.metadata.tonly$REF, indels.metadata.tonly$ALT,
                sep=":")
panel.str = paste(panel.nv.total$CHROM, panel.nv.total$POS,
                  panel.nv.total$REF, panel.nv.total$ALT,
                  sep=":")

MAX.THRESH = 5  # originally =10
snvs.in.panel.idx = vector(mode="list", length=MAX.THRESH)
indels.in.panel.idx = vector(mode="list", length=MAX.THRESH)

for (min.nv in 1:MAX.THRESH) {
    # Identify somatic variants colliding with variants present in the host panel
    idx = panel.nv.total$NV_TOTAL >= min.nv
    snvs.in.panel.idx[[min.nv]] = snv.str %in% panel.str[idx]
    indels.in.panel.idx[[min.nv]] = ind.str %in% panel.str[idx]
}

trinucs = substr(str_split_fixed(snvs.metadata.tonly$INFO, ";", 20)[, 12], 13, 15)
trinuc.N.idx = grepl("N", trinucs)  # SNVs with an 'N' in their context will be removed

# # For each threshold of supporting reads across the host panel (1-5), plot spectra
# cat("\nPlotting spectra of colliding variant sets...\n")
# # Make variant table for sigfit
# nrows = sum(sapply(snvs.in.panel.idx, function(idx) {
#     sum(idx & !trinuc.N.idx)
# }))
# var.table = matrix(NA, nrow=nrows, ncol=4)
# i = 1
# for (min.nv in 1:MAX.THRESH) {
#     idx = snvs.in.panel.idx[[min.nv]] & !trinuc.N.idx
#     n = sum(idx)
#     var.table[i:(i+n-1), ] = cbind(paste("Variants with", min.nv, "or more supporting reads across host panel"),
#                              snvs.metadata.tonly$REF[idx],
#                              snvs.metadata.tonly$ALT[idx],
#                              trinucs[idx])
#     i = i + n
# }
# # Plot spectra
# var.spectra = build_catalogues(var.table)
# invisible(plot_spectrum(var.spectra, pdf_path=paste0("../output/", Sys.Date(), "_Spectra_HostPanelVariants.pdf")))
# cat("Spectra plotted to output/ folder\n")

# Select variants with <5 total supporting reads across the host panel
MIN.NV = 5
not.in.panel.snv.idx = !snvs.in.panel.idx[[MIN.NV]] & !trinuc.N.idx
not.in.panel.ind.idx = !indels.in.panel.idx[[MIN.NV]]


# 5) Repeat filter
# Discard variants in, or within 5 bp of, repeat regions (simple/low-complex./tandem repeats, excluding exons)
# Load repeat coordinates
# # Processing repeat regions: union of RepeatMasker and TRF, subtraction of Ensembl exons (done in local)
# cf3.rmsk = read.table("CanFam3.1_RepeatMasker_UCSC.tsv", header=T, comment.char="", stringsAsFactors=F)
# cf3.rmsk = cf3.rmsk[cf3.rmsk$genoName %in% paste0("chr", c(1:38, "X")) & 
#                         cf3.rmsk$repClass %in% c("Low_complexity", "Simple_repeat"), ]
# cf3.rmsk$genoName = substring(cf3.rmsk$genoName, 4)
# cf3.rmsk.gr = makeGRangesFromDataFrame(cf3.rmsk, seqnames="genoName", start="genoStart", end="genoEnd",
#                                        ignore.strand=T, starts.in.df.are.0based=T)
# cf3.trf = read.table("CanFam3.1_SimpleRepeat_UCSC.tsv", header=T, comment.char="", stringsAsFactors=F)
# cf3.trf = cf3.trf[cf3.trf$chrom %in% paste0("chr", c(1:38, "X")), ]
# cf3.trf$chrom = substring(cf3.trf$chrom, 4)
# cf3.trf.gr = makeGRangesFromDataFrame(cf3.trf, seqnames="chrom", start="chromStart", end="chromEnd",
#                                       ignore.strand=T, starts.in.df.are.0based=T)
# cf3.repeats.gr = GenomicRanges::union(cf3.rmsk.gr, cf3.trf.gr)
# load("~/Desktop/R/DAT_2017-11-29_EnsemblAnnotation_Genes-Exons-UTRs.RData")  # 1-based start
# cf3.repeats.minus.exons.gr = GenomicRanges::setdiff(cf3.repeats.gr, exons.all.gr)
# save(cf3.repeats.gr, cf3.repeats.minus.exons.gr, file="DAT_2018-11-08_CanFam3.1_Repeats_Simple-LowCompl-TRF_UCSC.RData")
load("DAT_2018-11-08_CanFam3.1_Repeats_Simple-LowCompl-TRF_UCSC.RData")

# Include 5-bp buffer around region coordinates
BUFFER = 5
start(cf3.repeats.minus.exons.gr) = start(cf3.repeats.minus.exons.gr) - BUFFER
end(cf3.repeats.minus.exons.gr) = end(cf3.repeats.minus.exons.gr) + BUFFER

# Select variants that do not overlap repeats
snvs.metadata.tonly.gr = makeGRangesFromDataFrame(snvs.metadata.tonly, 
                                                  seqnames="CHROM", start.field="POS", end="POS", ignore.strand=T)
indel.len = sapply(indels.metadata.tonly$REF, function(ref) max(nchar(ref) - 1, 1))
indels.metadata.tonly.gr = makeGRangesFromDataFrame(indels.metadata.tonly, 
                                                    seqnames="CHROM", start="POS", end="POS", ignore.strand=T,
                                                    starts.in.df.are.0based=T)  # to account for inclusion of base before indel
end(indels.metadata.tonly.gr) = end(indels.metadata.tonly.gr) + indel.len

not.in.repeat.snv.idx = is.na(findOverlaps(snvs.metadata.tonly.gr, cf3.repeats.minus.exons.gr, select="first"))
not.in.repeat.ind.idx = is.na(findOverlaps(indels.metadata.tonly.gr, cf3.repeats.minus.exons.gr, select="first"))


# 6) Strand bias filter
# Discard variants with ≤20% of total supporting reads on either strand
# (It follows from filter (1) that the minimum total number of supporting reads is NF+NR=5)
nf.nr.snv = apply(str_split_fixed(snvs.metadata.tonly$INFO, ";", 20)[, 8:9], 2,
                  function(x) as.integer(str_split_fixed(x, "=", 2)[, 2]))
nf.nr.ind = apply(str_split_fixed(indels.metadata.tonly$INFO, ";", 20)[, 8:9], 2,
                  function(x) as.integer(str_split_fixed(x, "=", 2)[, 2]))

# Select variants with >20% of reads on each strand
MIN.RATIO = 0.2
no.bias.snv.idx = apply(nf.nr.snv, 1, function(nf.nr) {
    min(nf.nr) / sum(nf.nr) > MIN.RATIO
})
no.bias.ind.idx = apply(nf.nr.ind, 1, function(nf.nr) {
    min(nf.nr) / sum(nf.nr) > MIN.RATIO
})


# 7) Unmatched tumour-unique filter
# Discard variants which are unique to an unmatched tumour, and have VAF < 0.3 * purity
# NOTE: This filter can only be applied after obtaining purity estimates in Step 7,
#       and should be commented when running this script for the first time.

# Load tumour purity estimates
load("TumourPurity.RData")

# Match tumours and hosts by number to identify unmatched tumours
hostmatch <- fread(commandArgs(TRUE)[1], na.strings = c("", "NA", "N/A"))
setnames(hostmatch, c("Tumour", "Host"))
hostmatch[data.table(Tumour=samples[tumours]), , on = "Tumour"]
unmatched <- hostmatch[data.table(Tumour=samples[tumours]), is.na(Host), on = "Tumour"]

# Identify tumour-unique variants with VAF < 0.3*purity in each unmatched tumour
MIN.VAF = 0.3
unmatched.low.vaf.snv.idx = sapply(which(unmatched), function(i) {
    tumour.unique.idx[, i] &
       snvs.vaf.tonly[, tumours][, i] < MIN.VAF * purity[i]
})
unmatched.low.vaf.ind.idx = sapply(which(unmatched), function(i) {
    tumour.unique.indels.idx[, i] &
        indels.vaf.tonly[, tumours][, i] < MIN.VAF * purity[i]
})

not.unmatched.low.vaf.snv.idx = rowSums(unmatched.low.vaf.snv.idx) == 0
not.unmatched.low.vaf.ind.idx = rowSums(unmatched.low.vaf.ind.idx) == 0


# Apply filters
cat("Variants exclusively discarded by each filter:\n")
cat("Low support filter:", 
    num(sum(not.unmatched.low.vaf.snv.idx & no.bias.snv.idx & not.in.repeat.snv.idx & 
                not.in.panel.snv.idx & not.mismap.snv.idx & no.host.supp.snv.idx & !five.reads.snv.idx)), "SNVs and", 
    num(sum(not.unmatched.low.vaf.ind.idx & no.bias.ind.idx & not.in.repeat.ind.idx & 
                not.in.panel.ind.idx & not.mismap.ind.idx & no.host.supp.ind.idx & !five.reads.ind.idx)), "indels\n")
cat("Host support filter:", 
    num(sum(not.unmatched.low.vaf.snv.idx & no.bias.snv.idx & not.in.repeat.snv.idx & 
                not.in.panel.snv.idx & not.mismap.snv.idx & !no.host.supp.snv.idx & five.reads.snv.idx)), "SNVs and", 
    num(sum(not.unmatched.low.vaf.ind.idx & no.bias.ind.idx & not.in.repeat.ind.idx & 
                not.in.panel.ind.idx & not.mismap.ind.idx & !no.host.supp.ind.idx & five.reads.ind.idx)), "indels\n")
cat("Mismapping filter:", 
    num(sum(not.unmatched.low.vaf.snv.idx & no.bias.snv.idx & not.in.repeat.snv.idx & 
                not.in.panel.snv.idx & !not.mismap.snv.idx & no.host.supp.snv.idx & five.reads.snv.idx)), "SNVs and", 
    num(sum(not.unmatched.low.vaf.ind.idx & no.bias.ind.idx & not.in.repeat.ind.idx & 
                not.in.panel.ind.idx & !not.mismap.ind.idx & no.host.supp.ind.idx & five.reads.ind.idx)), "indels\n")
cat("Host panel filter:", 
    num(sum(not.unmatched.low.vaf.snv.idx & no.bias.snv.idx & not.in.repeat.snv.idx & 
                !not.in.panel.snv.idx & not.mismap.snv.idx & no.host.supp.snv.idx & five.reads.snv.idx)), "SNVs and", 
    num(sum(not.unmatched.low.vaf.ind.idx & no.bias.ind.idx & not.in.repeat.ind.idx & 
                !not.in.panel.ind.idx & not.mismap.ind.idx & no.host.supp.ind.idx & five.reads.ind.idx)), "indels\n")
cat("Simple/low-complexity/tandem repeat filter:", 
    num(sum(not.unmatched.low.vaf.snv.idx & no.bias.snv.idx & !not.in.repeat.snv.idx & 
                not.in.panel.snv.idx & not.mismap.snv.idx & no.host.supp.snv.idx & five.reads.snv.idx)), "SNVs and", 
    num(sum(not.unmatched.low.vaf.ind.idx & no.bias.ind.idx & !not.in.repeat.ind.idx & 
                not.in.panel.ind.idx & not.mismap.ind.idx & no.host.supp.ind.idx & five.reads.ind.idx)), "indels\n")
cat("Strand bias filter:", 
    num(sum(not.unmatched.low.vaf.snv.idx & !no.bias.snv.idx & not.in.repeat.snv.idx & 
                not.in.panel.snv.idx & not.mismap.snv.idx & no.host.supp.snv.idx & five.reads.snv.idx)), "SNVs and", 
    num(sum(not.unmatched.low.vaf.ind.idx & !no.bias.ind.idx & not.in.repeat.ind.idx & 
                not.in.panel.ind.idx & not.mismap.ind.idx & no.host.supp.ind.idx & five.reads.ind.idx)), "indels\n")
cat("Unmatched tumour-unique filter:", 
    num(sum(!not.unmatched.low.vaf.snv.idx & no.bias.snv.idx & not.in.repeat.snv.idx & 
                not.in.panel.snv.idx & not.mismap.snv.idx & no.host.supp.snv.idx & five.reads.snv.idx)), "SNVs and", 
    num(sum(!not.unmatched.low.vaf.ind.idx & no.bias.ind.idx & not.in.repeat.ind.idx & 
                not.in.panel.ind.idx & not.mismap.ind.idx & no.host.supp.ind.idx & five.reads.ind.idx)), "indels\n")

cat("\nApplying filters...\n")
filter.snv.idx = not.unmatched.low.vaf.snv.idx & no.bias.snv.idx & not.in.repeat.snv.idx & 
    not.in.panel.snv.idx & not.mismap.snv.idx & no.host.supp.snv.idx & five.reads.snv.idx
snvs.metadata.tonly = snvs.metadata.tonly[filter.snv.idx, ]
snvs.nr.tonly = snvs.nr.tonly[filter.snv.idx, ]
snvs.nv.tonly = snvs.nv.tonly[filter.snv.idx, ]
snvs.vaf.tonly = snvs.vaf.tonly[filter.snv.idx, ]
tumour.unique.idx = tumour.unique.idx[filter.snv.idx, ]
tumour.unique.any.idx = tumour.unique.any.idx[filter.snv.idx]

filter.ind.idx = not.unmatched.low.vaf.ind.idx & no.bias.ind.idx & not.in.repeat.ind.idx & 
    not.in.panel.ind.idx & not.mismap.ind.idx & no.host.supp.ind.idx & five.reads.ind.idx
indels.metadata.tonly = indels.metadata.tonly[filter.ind.idx, ]
indels.nr.tonly = indels.nr.tonly[filter.ind.idx, ]
indels.nv.tonly = indels.nv.tonly[filter.ind.idx, ]
indels.vaf.tonly = indels.vaf.tonly[filter.ind.idx, ]
tumour.unique.indels.idx = tumour.unique.indels.idx[filter.ind.idx, ]
tumour.unique.indels.any.idx = tumour.unique.indels.any.idx[filter.ind.idx]

cat("After filtering:", num(nrow(snvs.metadata.tonly)), "tumour-only SNVs and",
    num(nrow(indels.metadata.tonly)), "tumour-only indels\n")


# Save objects
cat("Saving indices and tables to data/ directory...\n")
dir.create("data")
save(samples, tumours, hosts, excluded, unmatched, MIN.READS, tumour.unique.idx, 
     tumour.unique.any.idx, tumour.unique.indels.idx, tumour.unique.indels.any.idx, 
     file="data/VariantTables_Indices_Filt.RData")
saveRDS(snvs.metadata.tonly, file="data/VariantTables_Split_Filt_snvs_Metadata.tonly.RDS")
saveRDS(snvs.nr.tonly, file="data/VariantTables_Split_Filt_snvs_NR.tonly.RDS")
saveRDS(snvs.nv.tonly, file="data/VariantTables_Split_Filt_snvs_NV.tonly.RDS")
saveRDS(snvs.vaf.tonly, file="data/VariantTables_Split_Filt_snvs_VAF.tonly.RDS")
saveRDS(indels.metadata.tonly, file="data/VariantTables_Split_Filt_indels_Metadata.tonly.RDS")
saveRDS(indels.nr.tonly, file="data/VariantTables_Split_Filt_indels_NR.tonly.RDS")
saveRDS(indels.nv.tonly, file="data/VariantTables_Split_Filt_indels_NV.tonly.RDS")
saveRDS(indels.vaf.tonly, file="data/VariantTables_Split_Filt_indels_VAF.tonly.RDS")
