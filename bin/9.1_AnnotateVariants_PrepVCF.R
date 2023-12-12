#!/usr/bin/env Rscript
# Analysis of CTVT whole-genome variants
# Step 9. Output variant sets in VCF format for annotation using VEP
# Adrian Baez-Ortega, 2018

# Load variant data
snvs.metadata.host    <- readRDS("VariantTables_Split_snvs_Metadata.host.RDS")
indels.metadata.host  <- readRDS("VariantTables_Split_indels_Metadata.host.RDS")
snvs.metadata.tonly   <- readRDS("VariantTables_Split_Filt_snvs_Metadata.tonly.RDS")
indels.metadata.tonly <- readRDS("VariantTables_Split_Filt_indels_Metadata.tonly.RDS")


# Erase QUAL, FILTER and INFO fields and add sequential ID
snvs.metadata.host[, 6:8] = "."
snvs.metadata.tonly[, 6:8] = "."
indels.metadata.host[, 6:8] = "."
indels.metadata.tonly[, 6:8] = "."
snvs.metadata.host$ID = seq_len(nrow(snvs.metadata.host))
snvs.metadata.tonly$ID = seq_len(nrow(snvs.metadata.tonly))
indels.metadata.host$ID = seq_len(nrow(indels.metadata.host))
indels.metadata.tonly$ID = seq_len(nrow(indels.metadata.tonly))


# Output
dir.create("data")
write.table(snvs.metadata.host, file="data/HostSNVs_ToAnnotate.vcf", sep="\t", quote=F, row.names=F, col.names=F)
write.table(snvs.metadata.tonly, file="data/SomaticSNVs_ToAnnotate.vcf", sep="\t", quote=F, row.names=F, col.names=F)
write.table(indels.metadata.host, file="data/HostIndels_ToAnnotate.vcf", sep="\t", quote=F, row.names=F, col.names=F)
write.table(indels.metadata.tonly, file="data/SomaticIndels_ToAnnotate.vcf", sep="\t", quote=F, row.names=F, col.names=F)

cat("\nVCF files for annotation of SNVs and indels written in data/ directory.\n")
