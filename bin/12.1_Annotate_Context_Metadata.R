#!/usr/bin/env Rscript
# Examines the split variant tables and assigns somatic / germline and
# removed status

library(data.table)

# Load tables
snvs.metadata       <- readRDS("VariantTables_Filtered_snvs_Metadata.RDS")
snvs.metadata.host  <- readRDS("VariantTables_Split_snvs_Metadata.host.RDS")
snvs.metadata.tonly <- readRDS("VariantTables_Split_Filt_snvs_Metadata.tonly.RDS")

# Add an original Index to the metadata,
# so we never lose the original sorting order
# - essential as the separate tables snvs.nr
# and snvs.nv MUST stay in sync
snvs.metadata$original_index <- seq(1, nrow(snvs.metadata))

setDT(snvs.metadata)
setDT(snvs.metadata.host)
setDT(snvs.metadata.tonly)
snvs.metadata[,
              CONTEXT := stringr::str_match(INFO,
                                            "SC=\\w{9}(\\w{3})\\w{9};")[, 2]]

all_snvs <- snvs.metadata[, .(CHROM, POS, REF, ALT, CONTEXT, original_index)]
host <- snvs.metadata.host[, .(CHROM, POS, REF, ALT)]
filtered_somatic <- snvs.metadata.tonly[, .(CHROM, POS, REF, ALT)]
rm(snvs.metadata)
rm(snvs.metadata.host)
rm(snvs.metadata.tonly)

setkey(all_snvs, CHROM, POS, REF, ALT)
setkey(host)
setkey(filtered_somatic)

# Identify each row in 'all_snvs' as is_germline=TRUE if it also appears in
# 'host', and as is_somatic=TRUE if it appears in filtered_somatic
all_snvs[, c("is_germline", "is_somatic", "is_somatic_filter_fail") := FALSE]
all_snvs[host, is_germline := TRUE]
all_snvs[filtered_somatic, is_somatic := TRUE]

# Any remaining rows not in either category is a row that was at one point
# considered somatic, but later failed some filter.
all_snvs[is_germline == FALSE & is_somatic == FALSE,
         is_somatic_filter_fail := TRUE]

setorder(all_snvs, original_index)
all_snvs[, original_index := NULL]
snvs.metadata.annotated <- as.data.frame(all_snvs)
dir.create("data")
saveRDS(snvs.metadata.annotated,
     file = file.path("data",
        "VariantTables_Filtered_snvs_Annotated_Metadata.RDS"))

# Load tables
indels.metadata       <- readRDS("VariantTables_Filtered_indels_Metadata.RDS")
indels.metadata.host  <- readRDS("VariantTables_Split_indels_Metadata.host.RDS")
indels.metadata.tonly <- readRDS("VariantTables_Split_Filt_indels_Metadata.tonly.RDS")

# Add an original Index to the metadata,
# so we never lose the original sorting order
# - essential as the separate tables indels.nr
# and indels.nv MUST stay in sync
indels.metadata$original_index <- seq(1, nrow(indels.metadata))

setDT(indels.metadata)
setDT(indels.metadata.host)
setDT(indels.metadata.tonly)
indels.metadata[,
              CONTEXT := stringr::str_match(INFO,
                                            "SC=\\w{9}(\\w{3})\\w{9};")[, 2]]

all_indels <- indels.metadata[, .(CHROM, POS, REF, ALT, CONTEXT, original_index)]
host <- indels.metadata.host[, .(CHROM, POS, REF, ALT)]
filtered_somatic <- indels.metadata.tonly[, .(CHROM, POS, REF, ALT)]
rm(indels.metadata)
rm(indels.metadata.host)
rm(indels.metadata.tonly)

setkey(all_indels, CHROM, POS, REF, ALT)
setkey(host)
setkey(filtered_somatic)

# Identify each row in 'all_indels' as is_germline=TRUE if it also appears in
# 'host', and as is_somatic=TRUE if it appears in filtered_somatic
all_indels[, c("is_germline", "is_somatic", "is_somatic_filter_fail") := FALSE]
all_indels[host, is_germline := TRUE]
all_indels[filtered_somatic, is_somatic := TRUE]

# Any remaining rows not in either category is a row that was at one point
# considered somatic, but later failed some filter.
all_indels[is_germline == FALSE & is_somatic == FALSE,
         is_somatic_filter_fail := TRUE]

setorder(all_indels, original_index)
all_indels[, original_index := NULL]
indels.metadata.annotated <- as.data.frame(all_indels)
saveRDS(indels.metadata.annotated,
     file = file.path("data",
        "VariantTables_Filtered_indels_Annotated_Metadata.RDS"))
