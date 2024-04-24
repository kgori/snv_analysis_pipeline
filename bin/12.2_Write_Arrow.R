#!/usr/bin/env Rscript
#' Split the data into separate tables for each sample
suppressPackageStartupMessages({
    library(data.table)
    library(arrow)
})

hostmatch_file <- commandArgs(trailing=TRUE)[1]

cat("Loading and combining data\n")
snvs.metadata.annotated   <- readRDS("VariantTables_Filtered_snvs_Annotated_Metadata.RDS")
indels.metadata.annotated <- readRDS("VariantTables_Filtered_indels_Annotated_Metadata.RDS")
snvs.metadata.annotated$is_indel <- FALSE
indels.metadata.annotated$is_indel <- TRUE
variants.metadata.annotated <- rbind(snvs.metadata.annotated,
                                     indels.metadata.annotated)
rm(snvs.metadata.annotated, indels.metadata.annotated)

snvs.nr     <- readRDS("VariantTables_Filtered_snvs_NR.RDS")
indels.nr   <- readRDS("VariantTables_Filtered_indels_NR.RDS")
variants.nr <- rbind(snvs.nr, indels.nr)
rm(snvs.nr, indels.nr)

snvs.nv     <- readRDS("VariantTables_Filtered_snvs_NV.RDS")
indels.nv   <- readRDS("VariantTables_Filtered_indels_NV.RDS")
variants.nv <- rbind(snvs.nv, indels.nv)
rm(snvs.nv, indels.nv)

host_matches <- fread(hostmatch_file,
    na.strings = c("NA", "N/A", "na", "n/a", ""))
setnames(host_matches, c("tumour", "host"))
host_matches[host == "", host := NA]

setDT(variants.nr)
setDT(variants.nv)
setDT(variants.metadata.annotated)

outdir <- "data/dataset"
if (!dir.exists(outdir)) {
    dir.create(outdir, recursive = TRUE)
}

if (!dir.exists(outdir)) {
    stop("Could not create output directory")
}

# dir.create("data/fst")

tumours <- host_matches[, unique(tumour)]
if (!all(tumours %in% colnames(variants.nv))) {
    stop("Not all tumours are present in the data")
}

for (samplename_ in tumours) {
    cat(paste("Extracting sample", samplename_, "\n"))
    matched_host <- host_matches[tumour == samplename_, host]
    dt <- copy(variants.metadata.annotated)
    dt[, T_total_reads := variants.nr[, .SD, .SDcols = samplename_][[1]]]
    dt[, T_alt_reads := variants.nv[, .SD, .SDcols = samplename_][[1]]]
    dt[, T_vaf := T_alt_reads / T_total_reads]

    cat(paste("Adding annotations", "\n"))
    if (!is.na(matched_host)) {
        dt[, H_total_reads := variants.nr[, .SD, .SDcols = matched_host][[1]]]
        dt[, H_alt_reads := variants.nv[, .SD, .SDcols = matched_host][[1]]]
        dt[, H_vaf := H_alt_reads / H_total_reads]
    } else {
        dt[, H_total_reads := NA_integer_]
        dt[, H_alt_reads := NA_integer_]
        dt[, H_vaf := NA_real_]
    }
    dt[, samplename := samplename_]
    if (!is.na(matched_host)) {
        dt[, hostname := matched_host]
    } else {
        dt[, hostname := NA_character_]
    }

    setkey(dt, samplename, CHROM, POS, REF, ALT)

    cat(paste("Writing output to", outdir, "\n"))
    arrow::write_dataset(dt, outdir, partitioning = c("samplename", "CHROM"))
    rm(dt)
    gc()
    # write.fst(dt, paste0("data/fst/Sample_", samplename_, ".fst"))
}

