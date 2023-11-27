#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(logging))
basicConfig("INFO")

check_file <- function(filename) {
    if (!file.exists(filename)) {
        logerror(paste("File", filename, "not found"))
        stop("ERROR - file error")
    } else {
        return(normalizePath(filename))
    }
}

variant_depth_filter <- 
  function (dt, variant_column_names, min_variant_reads, min_vaf) {
    mask <- rep(FALSE, nrow(dt))
    if (length(min_variant_reads) == length(variant_column_names)) {
        thresholds <- data.frame(name = variant_column_names,
            threshold = min_variant_reads, stringsAsFactors = FALSE)
    }
    else {
        thresholds <- data.frame(name = variant_column_names,
            threshold = rep(min_variant_reads[1], length(variant_column_names)),
            stringsAsFactors = FALSE)
    }
    for (i in 1:nrow(thresholds)) {
        nv <- thresholds[i, ][["name"]]
        vaf <- sub(".nv", ".vaf", nv)
        min_nv <- thresholds[i, ][["threshold"]]
        mask <- mask | dt[, eval(as.name(nv)) >= min_nv | eval(as.name(vaf)) >=
            min_vaf]
    }
    mask[is.na(mask)] <- FALSE
    mask
}

# Compute VAFs
compute_vaf <- function(var, tot) {
    vaf <- var / tot
    vaf[is.na(vaf)] <- 0
    vaf
}

args <- commandArgs(trailing = TRUE)
if (length(args) < 1) {
    stop(paste("Usage: Rscript assign_germline.R <snvs_tsv>",
               "<output_tsv> [<optional:host_panel_snvs>]"))
}

snvs_tsv <- check_file(args[1])
smpl_assn_file <- check_file(args[2])
outfile <- args[3]

use_host_panel <- FALSE
if (length(args) == 4) {
    use_host_panel <- TRUE
    host_panel_file <- check_file(args[4])
    loginfo(sprintf("Using host panel %s\n", host_panel_file))
}

smpl_assn <- fread(smpl_assn_file)

loginfo(paste("Loading snvs from", snvs_tsv))
snvs <- fread(snvs_tsv, showProgress = TRUE)
colnames(snvs)[1] <- "CHROM"
snvs[, CHROM := as.character(CHROM)]
nr_cols <- grep("*\\.nr$", colnames(snvs), value = TRUE)

expected_host_nr <- paste(smpl_assn[Type == "H", Sample],
                          "nr", sep = ".")
host_nr <- intersect(expected_host_nr, nr_cols)
expected_tumour_nr <- paste(smpl_assn[Type == "T", Sample],
                            "nr", sep = ".")
tumour_nr <- intersect(expected_tumour_nr, nr_cols)

for (sample in c(host_nr, tumour_nr)) {
    loginfo(paste("Computing VAF for", sample))
    vaf_col <- sub(".nr", ".vaf", sample)
    variants_col <- sub(".nr", ".nv", sample)
    total_col <- sample
    snvs[, (vaf_col) := compute_vaf(eval(as.name(variants_col)),
                                    eval(as.name(total_col)))]
}

# Filter by germline panel
if (use_host_panel) {
    loginfo("Processing host panel file")

    # fread can read gzipped inputs directly
    panel <- fread(host_panel_file)

    # Find the overlap (SNVs are in data.table "snvs", provided elsewhere)
    panel[, I := 1] # indicator for later filtering
    setkey(panel, CHROM, POS, REF, ALT)

    panel_mask <- panel[snvs[, 1:4], !is.na(I)]
    rm(panel)
} else {
    panel_mask <- rep(FALSE, nrow(snvs))
}

loginfo("Identifying germline variants")
host_nv <- sub(".nr", ".nv", host_nr)
tumour_nv <- sub(".nr", ".nv", tumour_nr)
germline_mask <- variant_depth_filter(snvs, host_nv, 3, 0.05)
somatic_mask <- variant_depth_filter(snvs, tumour_nv, 5, 0.05)

snvs[, type := "unknown"]
snvs[panel_mask | germline_mask, type := "germline"]
snvs[somatic_mask & !panel_mask & !germline_mask, type := "somatic"]
snvs[is.na(type), type := "unknown"]

loginfo(paste("Writing results to", outfile))
fwrite(snvs, file = outfile, sep = '\t', showProgress = TRUE)

# Purity estimates
plotfile <- paste0(outfile, ".pdf")
loginfo(paste("Plotting purity estimates to", plotfile))
tumour_vafs <- sub(".nr", ".vaf", tumour_nr)
pdf(plotfile, onefile = TRUE)

purity <- rep(0, length(tumour_vafs))
names(purity) <- tumour_vafs

for (tumour_vaf in tumour_vafs) {
    loginfo(paste("Plotting purity of sample", tumour_vaf))
    vafs <- snvs[type == "somatic" & eval(as.name(sub(".vaf",".nv", tumour_vaf))) >= 5, eval(as.name(tumour_vaf))]
    vafs_centred <- vafs[vafs > 0.0 & vafs < 1]
    peaks <- modes::amps(vafs_centred)$Peaks
    het_peak <- peaks[which.max(peaks[,2]), 1]
    purity_estimate <- min(1, 2 * het_peak)
    purity[[tumour_vaf]] <- purity_estimate
    hist(vafs, breaks = seq(0, 1, length.out = 51), col = "gray54", border = "white",
         main = paste0(tumour_vaf, " purity = ", 100*round(purity_estimate, 3), "%"),
         xlab = "VAF (somatic variants)")
    abline(v = het_peak, col = "firebrick2", lwd = 2.5)
}
dev.off()

loginfo("Saving purity estimates to TumourPurity.RData")
save(purity, file = "TumourPurity.RData")

