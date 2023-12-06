#!/usr/bin/env Rscript
# Analysis of CTVT whole-genome variants
# Step 5. Output somatic variants VCF for genotyping on the low-coverage host panel
# Adrian Baez-Ortega, 2018

# Load variant data
snvs.metadata.tonly <- readRDS("VariantTables_Split_snvs_Metadata.tonly.RDS")
indels.metadata.tonly <- readRDS("VariantTables_Split_indels_Metadata.tonly.RDS")

all.metadata.tonly = rbind(snvs.metadata.tonly, indels.metadata.tonly)
all.metadata.tonly[, 6:8] = "."

if (!dir.exists("data")) {
  dir.create("data")
}

write.table(all.metadata.tonly,
    file = "data/SomaticVariants_ToGenotype.vcf",
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = FALSE)

cat("\nVCF file of somatic SNVs and indels written in data/ directory.\n")
