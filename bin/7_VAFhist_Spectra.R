#!/usr/bin/env Rscript
# Analysis of CTVT whole-genome variants
# Step 7. Exploratory variant analyses: coverage and VAF histograms, mutational spectra
# Adrian Baez-Ortega, 2018

# Load packages
library(data.table)
library(stringr)
library(sigfit, quietly=T)

# Load variant data
load("VariantTables_Indices_Filt.RData")
snvs.metadata.host <-  readRDS("VariantTables_Split_snvs_Metadata.host.RDS")
snvs.nr.host <-        readRDS("VariantTables_Split_snvs_NR.host.RDS")
snvs.nv.host <-        readRDS("VariantTables_Split_snvs_NV.host.RDS")
snvs.vaf.host <-       readRDS("VariantTables_Split_snvs_VAF.host.RDS")
snvs.metadata.tonly <- readRDS("VariantTables_Split_Filt_snvs_Metadata.tonly.RDS")
snvs.nr.tonly <-       readRDS("VariantTables_Split_Filt_snvs_NR.tonly.RDS")
snvs.nv.tonly <-       readRDS("VariantTables_Split_Filt_snvs_NV.tonly.RDS")
snvs.vaf.tonly <-      readRDS("VariantTables_Split_Filt_snvs_VAF.tonly.RDS")

# Load CanFam3.1 whole-genome trinucleotide frequencies (from Biostrings::trinucleotideFrequency)
load("DAT_2018-11-12_CanFam3.1_TrinucFreqs_WG.RData")

# Reverse complement function (for string vectors)
rev.comp = function(nucleotide.list) {
    sapply(nucleotide.list, function(nucleotides) {
        paste(
            rev(sapply(strsplit(nucleotides, "")[[1]], function(nuc) {
                if (nuc == "A") "T"
                else if (nuc == "C") "G"
                else if (nuc == "G") "C"
                else if (nuc == "T") "A"
            })),
            collapse="")
    }, USE.NAMES=F)
}

# Output file paths
dir.create("output")
dir.create("data")
OUTPUT = list(
    COV.HIST.T = paste0("output/", Sys.Date(), "_CoverageHistograms_Tumours.pdf"),
    COV.HIST.H = paste0("output/", Sys.Date(), "_CoverageHistograms_Hosts.pdf"),
    VAF.HIST.T = paste0("output/", Sys.Date(), "_VAFhistograms_Tumours.pdf"),
    VAF.HIST.H = paste0("output/", Sys.Date(), "_VAFhistograms_Hosts.pdf"),
    VAR.SUM.T = paste0("output/", Sys.Date(), "_VariantCategories_Summary_Tumours.tsv"),
    SPECTRA.UNNORM = paste0("output/", Sys.Date(), "_Spectra_VariantGroups_Unnormalised.pdf"),
    SPECTRA.NORM = paste0("output/", Sys.Date(), "_Spectra_VariantGroups_Normalised.pdf"),
    SPECTRA = "data/VariantSpectra_Filt.RData",
    PURITY = "data/TumourPurity.RData")
    # QC.SUM.H = "../output/QC_Summary_Hosts.tsv"
    # QC.SUM.T = "../output/QC_Summary_Tumours.tsv"


# Obtain tumour purity estimates (from T-only VAF histogram modes)
# For each tumour, obtain the VAF value corresponding to the mode of 
# the tumour VAF histogram for present T-only variants
cat("Calculating tumour purity estimates...\n")
purity = sapply(which(tumours), function(i) {
    present = snvs.nv.tonly[, i] >= MIN.READS
    hist.data = hist(snvs.vaf.tonly[present, i], breaks=seq(0, 1, 0.01), plot=F)
    2 * hist.data$mids[which.max(hist.data$counts)]
})
purity[purity > 1.0] <- 1.0
names(purity) = samples[tumours]


# Function to plot VAF/coverage histograms
plot.histogram = function(filename, index, nv.table, plot.table, max.x, keyword, xlab) {  # max.y
    cairo_pdf(filename, 10, 8, onefile=T)
    par(mar=c(5.1, 4.8, 4.1, 0.8))
    
    for (i in which(index)) {
        present = nv.table[, i] >= MIN.READS
        h = hist(plot.table[present, i][plot.table[present, i] <= max.x],
                 breaks=seq(0, max.x, len=100),
                 col="dodgerblue4", border="white", 
                 xlim=c(0, max.x), xlab=xlab, ylab="SNVs",  # ylim=c(0, max.y)
                 main=paste0(keyword, " histogram for ", samples[i], "\n",
                             prettyNum(sum(present), big.mark=","), " SNVs (≥", MIN.READS, " supporting reads)"))
        if (keyword == "VAF") {
            #segments(x0=0.253, y0=0.13, y1=max.y, col="darkgrey", lwd=1.5)
            #text(x=0.106, y=0.98 * max.y, labels="Tumour contamination", col="grey40")
            abline(v=0.253, col="darkgrey", lwd=1.5)
            mtext("Tumour contamination", side=3, line=-6, at=0.106, col="grey40")
        }
    }
    invisible(dev.off())
}


# Plot coverage histograms per host and per tumour (only SNVs)
cat("Plotting coverage and VAF histograms...\n")
MAX.X = 300
# MAX.Y = 7000
plot.histogram(OUTPUT$COV.HIST.H, hosts, snvs.nv.host, snvs.nr.host, MAX.X, "Coverage", "Total read coverage")
plot.histogram(OUTPUT$COV.HIST.T, tumours, snvs.nv.tonly, snvs.nr.tonly, MAX.X, "Coverage", "Total read coverage")


# Plot VAF histograms of germline variants per host (only SNVs)
# MAX.Y = 50000
plot.histogram(OUTPUT$VAF.HIST.H, hosts, snvs.nv.host, snvs.vaf.host, 1, "VAF", "VAF in host")


# Plot VAF histograms of categorised variants per tumour (only SNVs)
# For each tumour, plot (on the same page) the VAF of variants that are:
#  - Tumour-only
#  - Tumour-unique
#  - Not in matched host (i.e. variant does not fit any of the categories above)
#  - Homozygous in matched host (VAF>0.75 in matched host)
#  - Heterozygous in matched host (0.25<VAF<0.75 in matched host)
#  - Contamination in matched host (VAF<0.25 in matched host)

# Identify tumours without matched host
hostmatch <- fread(commandArgs(TRUE)[1], na.strings = c("", "NA", "N/A"))
setnames(hostmatch, c("Tumour", "Host"))
hostmatch[data.table(Tumour=samples[tumours]), , on = "Tumour"]
unmatched <- hostmatch[data.table(Tumour=samples[tumours]), is.na(Host), on = "Tumour"]

# Plot histograms and collect variant numbers
cairo_pdf(OUTPUT$VAF.HIST.T, 14, 8, onefile=T)
par(mfrow=c(2, 3), mar=c(5.2, 2.1, 4, 0.5), oma=c(0, 3, 3, 0.5))
#MAX.Y = 50000

variant.categ = t(sapply(which(tumours), function(i) {
    
    # Get index of matched host
    tumour <- samples[i]
    matched_host_name <- hostmatch[Tumour == tumour, Host]
    if (is.na(matched_host_name)) {
        matched.host <- NA
    } else {
        matched.host <- which(samples == matched_host_name)
    }
    print(paste(tumour, matched_host_name))
    
    # Get present variants
    present.tonly = snvs.nv.tonly[, i] >= MIN.READS
    present.host = snvs.nv.host[, i] >= MIN.READS
    present.tunique = tumour.unique.idx[, samples[i]]
    
    # Plot variants in each category:
    # T-only
    hist(snvs.vaf.tonly[present.tonly, i], 
         breaks=seq(0, 1, 0.01), col="dodgerblue4", border="white", xpd=NA,
         xlab="VAF in tumour", ylab="SNVs", xlim=c(0, 1),  # ylim=c(0, MAX.Y)
         main=paste0(prettyNum(sum(present.tonly), big.mark=","), " tumour-only SNVs (≥", 
                     MIN.READS, " supporting reads)"))
    # T-unique
    hist(snvs.vaf.tonly[present.tunique, i], 
         breaks=seq(0, 1, 0.01), col="dodgerblue4", border="white",
         xlab="VAF in tumour", ylab="", xlim=c(0, 1),  # ylim=c(0, MAX.Y)
         main=paste0(prettyNum(sum(present.tunique), big.mark=","), " tumour-unique SNVs (≥", 
                     MIN.READS, " supporting reads)"))
    
    # If no matched host: plot nothing
    if (is.na(matched.host)) {
        host.avail = "N"
        host.indices = matrix(NA, nrow=1, ncol=4)
        for (j in 1:4) {
            hist(1, border=NA, xlab="VAF in tumour", ylab="", 
                 xlim=c(0, 1), main="(No matched host)")   # ylim=c(0, MAX.Y)
        }
    }
    
    # If matched host available: plot non-tumour-only variants
    else {
        host.avail = "Y"
        host.indices = cbind(
            # Found in tumour & not in matched host
            "not present" = present.host & snvs.vaf.host[, matched.host] == 0,
            
            # Found in tumour & homozygous in matched host
            "homozygous" = present.host & snvs.vaf.host[, matched.host] > 0.75,
            
            # Found in tumour & heterozygous in matched host
            "heterozygous" = present.host & snvs.vaf.host[, matched.host] > 0.25 & 
                snvs.vaf.host[, matched.host] < 0.75,
            
            # Found in tumour & contamination in matched host
            "contamination" = present.host & snvs.vaf.host[, matched.host] > 0 & 
                snvs.vaf.host[, matched.host] < 0.25
        )
        
        for (j in 1:ncol(host.indices)) {
            hist(snvs.vaf.host[host.indices[, j], i], 
                 breaks=seq(0, 1, 0.01), col="dodgerblue4", border="white", xpd=NA,
                 xlab="VAF in tumour", ylab=ifelse(j == 2, "SNVs", ""), xlim=c(0, 1),  # ylim=c(0, MAX.Y), 
                 main=paste0(prettyNum(sum(host.indices[, j]), big.mark=","), " SNVs (≥", 
                             MIN.READS, " supporting reads) that are\n", 
                             colnames(host.indices)[j], " in ", samples[matched.host]))
        }
    }
    title(paste("\nVAF histograms for", samples[i]), outer=T, cex.main=1.3)
    
    c(samples[i], purity[samples[i]], sum(present.tonly), 
      sum(present.tunique), colSums(host.indices), host.avail)
}))

invisible(dev.off())


# Output table of variants per tumour and category
colnames(variant.categ) = c("Sample", "Purity", "Tumour-only SNVs", "Tumour-unique SNVs",
                            "SNVs present in T and not in H",
                            "SNVs present in T and homozygous in H", 
                            "SNVs present in T and heterozygous in H",
                            "SNVs present in T and contamination in H",
                            "Matched H available")
write.table(variant.categ, row.names=F, quote=F, sep="\t",
            file=OUTPUT$VAR.SUM.T)


# Obtain tumour 'QC index', defined as the product (coverage * purity)
# qc.index = purity * coverage[tumours]
# # Output QC summaries for hosts and tumours
# cat("Writing sample QC summaries to output directory...\n")
# hosts.sum = cbind(samples[hosts],
#                   t(sapply(which(hosts), function(i) {
#                       c(sum(snvs.nv[, i] >= MIN.READS),
#                         sum(indels.nv[, i] >= MIN.READS))
#                   })),
#                   coverage[hosts])
# colnames(hosts.sum) = c("Sample", "SNVs", "Indels", "31+ coverage fraction")
# write.table(hosts.sum, quote=F, row.names=F, sep="\t",
#             file=OUTPUT$QC.SUM.H)
# tumours.sum = cbind(samples[tumours],
#                     t(sapply(which(tumours), function(i) {
#                         c(sum(snvs.nv[tumour.only.idx, i] >= MIN.READS),
#                           sum(tumour.unique.idx[, samples[i]]),
#                           sum(indels.nv[tumour.only.indels.idx, i] >= MIN.READS),
#                           sum(tumour.unique.indels.idx[, samples[i]]))
#                     })),
#                     coverage[tumours], 
#                     purity,
#                     qc.index)
# colnames(tumours.sum) = c("Sample", "Tumour-only SNVs", "Tumour-unique SNVs", 
#                           "Tumour-only indels", "Tumour-unique indels",
#                           "31+ coverage fraction", "Purity", "QC index")
# write.table(tumours.sum, quote=F, row.names=F, sep="\t", file=OUTPUT$QC.SUM.T)


# Produce mutational spectra of:
cat("Producing mutational spectra...\n")
var.table = matrix(NA, nrow=1e8, ncol=4)

# Host variants
trinucs = substr(str_split_fixed(snvs.metadata.host$INFO, ";", 20)[, 12], 13, 15)
idx = !grepl("N", trinucs)
cat(sum(!idx), "host variants with ambiguous context (ignored)\n")
i = 1
n = sum(idx)
var.table[i:(i+n-1), ] = cbind("All host (germline) variants",
                               snvs.metadata.host$REF[idx],
                               snvs.metadata.host$ALT[idx],
                               trinucs[idx])
# All tumour-only variants
trinucs = substr(str_split_fixed(snvs.metadata.tonly$INFO, ";", 20)[, 12], 13, 15)
i = i + n
n = length(trinucs)
var.table[i:(i+n-1), ] = cbind("All tumour-only variants",
                               snvs.metadata.tonly$REF,
                               snvs.metadata.tonly$ALT,
                               trinucs)
# Pre-divergence tumour-only variants (≥37/40 tumours)
MIN.TUM = floor(37 / 40 * length(samples[tumours]))
idx = rowSums(snvs.nv.tonly[, tumours] >= MIN.READS) >= MIN.TUM
i = i + n
n = sum(idx)
var.table[i:(i+n-1), ] = cbind("Pre-divergence tumour-only variants",
                               snvs.metadata.tonly$REF[idx],
                               snvs.metadata.tonly$ALT[idx],
                               trinucs[idx])
# All tumour-unique variants
i = i + n
n = sum(tumour.unique.any.idx)
var.table[i:(i+n-1), ] = cbind("All tumour-unique variants",
                               snvs.metadata.tonly$REF[tumour.unique.any.idx],
                               snvs.metadata.tonly$ALT[tumour.unique.any.idx],
                               trinucs[tumour.unique.any.idx])
# Tumour-unique variants per tumour
for (j in 1:ncol(tumour.unique.idx)) {
    i = i + n
    n = sum(tumour.unique.idx[, j])
    if (n > 0) {
        var.table[i:(i+n-1), ] = cbind(paste("Tumour-unique variants in", colnames(tumour.unique.idx)[j]),
                                       snvs.metadata.tonly$REF[tumour.unique.idx[, j]],
                                       snvs.metadata.tonly$ALT[tumour.unique.idx[, j]],
                                       trinucs[tumour.unique.idx[, j]])
    }
}
var.table = var.table[1:(i+n-1), ]

# Collapse trinucleotide frequencies into pyrimidine middle bases
# and build normalising frequency vector for spectra
freqs.pyr = cf3.freqs
pyr.idx = substr(names(cf3.freqs), 2, 2) %in% c("C", "T")
for (i in which(pyr.idx)) {
    freqs.pyr[i] = freqs.pyr[i] + freqs.pyr[ rev.comp(names(freqs.pyr)[i]) ]
}
freqs.pyr = freqs.pyr[pyr.idx]
freqs.pyr.C = freqs.pyr[grep(".C.", names(freqs.pyr))]
freqs.pyr.T = freqs.pyr[grep(".T.", names(freqs.pyr))]
freqs.spectrum = c(rep(freqs.pyr.C, 3), rep(freqs.pyr.T, 3))

# Produce unnormalised and normalised spectra
print("Building catalogues...")
spectra.unnorm = build_catalogues(var.table)
print("Building spectra...")
spectra.norm = t(apply(spectra.unnorm, 1, function(x) {
    z = x / freqs.spectrum
    z / sum(z)
}))
rownames(spectra.norm) = paste(rownames(spectra.norm), "- Normalised")
rownames(spectra.unnorm) = paste(rownames(spectra.unnorm), "- Unnormalised")
spectra.norm = rbind("CanFam3.1 trinucleotide frequencies (whole genome)"=freqs.spectrum,
                     spectra.norm)

# Plot spectra
print("Plotting spectra...")
invisible(plot_spectrum(spectra.unnorm, pdf_path=OUTPUT$SPECTRA.UNNORM))
invisible(plot_spectrum(spectra.norm, pdf_path=OUTPUT$SPECTRA.NORM))


# Save purity and spectrum objects
print("Saving purity and spectrum data...")
save(purity, file=OUTPUT$PURITY)
save(spectra.unnorm, spectra.norm, file=OUTPUT$SPECTRA)

cat("\nTumour purity and spectrum data saved to data/ directory",
    "\nOutput plots and variant summary in output/ directory\n")
