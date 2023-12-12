#!/usr/bin/env Rscript
# Analysis of CTVT whole-genome variants
# Step 4. Define initial set of tumour-only variants
# Adrian Baez-Ortega, 2018

# Candidate somatic (tumour-only) variants are defined as those having ≥3 reads
# in ≥1 tumour, and VAF <0.25 in all hosts.
# Germline (host+tumour) variants are defined as those with VAF >0.25 in any host.
# All variants must be either somatic or germline.
# VAF histograms of tumours and hosts are produced to identify contaminated
# or mislabelled samples. Tumour-only (candidate somatic) variants will later be
# genotyped across the panel of low-coverage hosts in order to deplete the
# remaining germline component.

# Load packages
library(data.table)
library(stringr)
library(logging); basicConfig()
args <- commandArgs(trailingOnly = TRUE)
smpl_assn_file <- args[1]
hostmatch_file <- args[2]
purities_file <- args[3]

if (!file.exists(smpl_assn_file)) {
    logerror("Sample association file not found: %s", smpl_assn_file)
    stop("ERROR - Missing file input")
}

smpl_assn <- fread(smpl_assn_file)
setnames(smpl_assn, c("Sample", "Type"))
if (!(all(smpl_assn$Type %in% c("H", "T")))) {
    logerror("Something is wrong with the sample file. Types other than H|T were supplied in Column 2")
    stop("ERROR - Bad Input")
}

if (!file.exists(hostmatch_file)) {
    logerror("Host match file not found: %s", hostmatch_file)
    stop("ERROR - Missing file input")
}

hostmatch <- fread(hostmatch_file, na.strings = c("", "NA", "N/A"))
setnames(hostmatch, c("Tumour", "Host"))
setkey(hostmatch, "Tumour")
if (!(all(hostmatch[!is.na(Tumour), Tumour] %in% smpl_assn[Type == "T", Sample]))) {
    logerror("Something is wrong with the host match file. Tumour IDs not found in sample file")
    stop("ERROR - Bad Input")
}

if (!(all(hostmatch[!is.na(Host), Host] %in% smpl_assn[Type == "H", Sample]))) {
    logerror("Something is wrong with the host match file. Host IDs not found in sample file")
    stop("ERROR - Bad Input")
}

if (!file.exists(purities_file)) {
    logerror("Purities file not found: %s", purities_file)
    stop("ERROR - Missing file input")
}
load(purities_file)
names(purity) <- sub("\\.vaf$", "", names(purity))
purity <- data.table(Tumour = names(purity), Purity = unname(purity))
setkey(purity, Tumour)

# Function to print numbers with commas
num = function(x) prettyNum(x, big.mark=",")

# Load variant data
cat("Loading variant tables...\n")
snvs.metadata <-   readRDS("VariantTables_Filtered_snvs_Metadata.RDS")
snvs.nr <-         readRDS("VariantTables_Filtered_snvs_NR.RDS")
snvs.nv <-         readRDS("VariantTables_Filtered_snvs_NV.RDS")
snvs.vaf <-        readRDS("VariantTables_Filtered_snvs_VAF.RDS")
indels.metadata <- readRDS("VariantTables_Filtered_indels_Metadata.RDS")
indels.nr <-       readRDS("VariantTables_Filtered_indels_NR.RDS")
indels.nv <-       readRDS("VariantTables_Filtered_indels_NV.RDS")
indels.vaf <-      readRDS("VariantTables_Filtered_indels_VAF.RDS")

loginfo("Colnames snvs.nr - %s", colnames(snvs.nr))
loginfo("Colnames snvs.nv - %s", colnames(snvs.nv))
loginfo("Colnames indels.nr - %s", colnames(indels.nr))
loginfo("Colnames indels.nv - %s", colnames(indels.nv))

# Create tumour and host indices
tumour_samples <- smpl_assn[Type == "T", Sample]
host_samples <- smpl_assn[Type == "H", Sample]
tumours = colnames(snvs.nr) %in% tumour_samples
hosts = colnames(snvs.nr) %in% host_samples
if (!all(tumours | hosts)) {
    logerror("Some samples are neither tumours nor hosts!")
    stop("ERROR - Unexpected Values")
}

if (!all((colnames(snvs.nv) %in% tumour_samples) == tumours)) {
    logerror("Column names do not match between SNV input tables")
    stop("ERROR - Bad Input")
}

if (!all((colnames(indels.nr) %in% tumour_samples) == tumours)) {
    logerror("Column names do not match between SNV and Indel input tables")
    stop("ERROR - Bad Input")
}

if (!all((colnames(indels.nv) %in% tumour_samples) == tumours)) {
    logerror("Column names do not match between Indel input tables")
    stop("ERROR - Bad Input")
}

cat(sum(tumours), "tumours,", sum(hosts), "hosts\n\n")


# Identify tumour-only variants
# Define 'special' tumours whose matched hosts have different IDs.
# We match 982T to host 1322_399H because 982T and 399T are the same tumour,
# and both 982T and 1322_399H were sequenced on X10
# (Project 1322 was sequenced on X10, while Project 1382 was sequenced on HiSeq).
samples <- colnames(snvs.nr)

# NB. For this sample set, we consider unmatched tumours in the identification
# of tumour-only variants, since they seem to be pure enough based on tumour VAF.
low_purity <- purity[Purity < 0.25, Tumour]

# Unmatched tumour samples lacking a matched host sample are not considered
# in the identification of tumour-only variants, as contaminating germline
# variants from the host dogs cannot be identified in these samples.
# Match tumours and hosts by number to identify unmatched tumours.
unmatched <- hostmatch[purity][is.na(Host) & Purity < 0.95, Tumour]

excluded <- (samples %in% low_purity | samples %in% unmatched)
tumours.except.excluded = tumours & !excluded

cat(sum(excluded), "tumour samples will be excluded from tumour-only variant identification:\n")
cat(samples[excluded], sep=", ")
cat("\n\n")

# Define indices for tumour-only SNVs and indels
# (variants with VAF <0.25 in all hosts and ≥3 reads in ≥1 tumour)
MIN.READS = 3
MIN.VAF.H = 0.25
tumour.only.idx = rowSums(snvs.vaf[, hosts] >= MIN.VAF.H) == 0 &
    rowSums(snvs.nv[, tumours.except.excluded] >= MIN.READS) > 0

tumour.only.indels.idx = rowSums(indels.vaf[, hosts] >= MIN.VAF.H) == 0 &
    rowSums(indels.nv[, tumours.except.excluded] >= MIN.READS) > 0

cat(num(sum(tumour.only.idx)), "tumour-only SNVs and",
    num(sum(tumour.only.indels.idx)), "tumour-only indels\n")


# Identify somatic SNVs/indels present in just one tumour ('tumour-unique')
tumour.unique.idx = snvs.nv[tumour.only.idx, tumours] >= MIN.READS
tumour.unique.idx[rowSums(tumour.unique.idx) > 1, ] = FALSE
tumour.unique.any.idx = as.logical(rowSums(tumour.unique.idx))

tumour.unique.indels.idx = indels.nv[tumour.only.indels.idx, tumours] >= MIN.READS
tumour.unique.indels.idx[rowSums(tumour.unique.indels.idx) > 1, ] = FALSE
tumour.unique.indels.any.idx = as.logical(rowSums(tumour.unique.indels.idx))

cat(num(sum(tumour.unique.any.idx)), "tumour-unique SNVs and",
    num(sum(tumour.unique.indels.any.idx)), "tumour-unique indels.\n")
cat("\nTumour-unique SNVs per tumour:\n")
for (j in 1:ncol(tumour.unique.idx)) {
    cat("   ", samples[tumours][j], ":  ", num(sum(tumour.unique.idx[,j])), " SNVs\n", sep="")
}


# Build separate tables for host and tumour-only variants
snvs.metadata.host = snvs.metadata[!tumour.only.idx, ]
snvs.nr.host = snvs.nr[!tumour.only.idx, ]
snvs.nv.host = snvs.nv[!tumour.only.idx, ]
snvs.vaf.host = snvs.vaf[!tumour.only.idx, ]

snvs.metadata.tonly = snvs.metadata[tumour.only.idx, ]
snvs.nr.tonly = snvs.nr[tumour.only.idx, ]
snvs.nv.tonly = snvs.nv[tumour.only.idx, ]
snvs.vaf.tonly = snvs.vaf[tumour.only.idx, ]

indels.metadata.host = indels.metadata[!tumour.only.indels.idx, ]
indels.nr.host = indels.nr[!tumour.only.indels.idx, ]
indels.nv.host = indels.nv[!tumour.only.indels.idx, ]
indels.vaf.host = indels.vaf[!tumour.only.indels.idx, ]

indels.metadata.tonly = indels.metadata[tumour.only.indels.idx, ]
indels.nr.tonly = indels.nr[tumour.only.indels.idx, ]
indels.nv.tonly = indels.nv[tumour.only.indels.idx, ]
indels.vaf.tonly = indels.vaf[tumour.only.indels.idx, ]


# Save objects
cat("Saving indices and tables to data/ directory...\n")
dir.create("data")
save(samples, tumours, hosts, excluded, MIN.READS, tumour.unique.idx,
     tumour.unique.any.idx, tumour.unique.indels.idx, tumour.unique.indels.any.idx,
     file="data/VariantTables_Indices.RData")
saveRDS(snvs.metadata.host,    file="data/VariantTables_Split_snvs_Metadata.host.RDS")
saveRDS(snvs.nr.host,          file="data/VariantTables_Split_snvs_NR.host.RDS")
saveRDS(snvs.nv.host,          file="data/VariantTables_Split_snvs_NV.host.RDS")
saveRDS(snvs.vaf.host,         file="data/VariantTables_Split_snvs_VAF.host.RDS")
saveRDS(indels.metadata.host,  file="data/VariantTables_Split_indels_Metadata.host.RDS")
saveRDS(indels.nr.host,        file="data/VariantTables_Split_indels_NR.host.RDS")
saveRDS(indels.nv.host,        file="data/VariantTables_Split_indels_NV.host.RDS")
saveRDS(indels.vaf.host,       file="data/VariantTables_Split_indels_VAF.host.RDS")
saveRDS(snvs.metadata.tonly,   file="data/VariantTables_Split_snvs_Metadata.tonly.RDS")
saveRDS(snvs.nr.tonly,         file="data/VariantTables_Split_snvs_NR.tonly.RDS")
saveRDS(snvs.nv.tonly,         file="data/VariantTables_Split_snvs_NV.tonly.RDS")
saveRDS(snvs.vaf.tonly,        file="data/VariantTables_Split_snvs_VAF.tonly.RDS")
saveRDS(indels.metadata.tonly, file="data/VariantTables_Split_indels_Metadata.tonly.RDS")
saveRDS(indels.nr.tonly,       file="data/VariantTables_Split_indels_NR.tonly.RDS")
saveRDS(indels.nv.tonly,       file="data/VariantTables_Split_indels_NV.tonly.RDS")
saveRDS(indels.vaf.tonly,      file="data/VariantTables_Split_indels_VAF.tonly.RDS")
