#!/usr/bin/env Rscript
# Analysis of CTVT whole-genome variants
# Step 3. Filter variants
# Adrian Baez-Ortega, 2018

# Filters:
# 1. Low-coverage filter: Exclude variants with median coverage across
#    hosts of <10 reads or >200 reads.
# 2. Host filter: Exclude variants that occur with 1–4 reads in ≥1 host,
#    and never occur with ≥5 reads in any sample (tumour or host).
# 3. Tumour filter: Exclude variants that occur with 1–4 reads in N≥2
#    tumours, and never occur with ≥5 reads in any tumour; except for
#    those variants that do not occur in all tumours (N<40) and for
#    which the N tumours with 1–4 reads are phylogenetically related
#    (i.e. the variant involves exactly 1 change in the tree).
# 4. Non-chromosomal filter: Exclude variants which are not in
#    chromosomes 1-38 or X.
# 5. Duplicated variant filter (we do not expect duplicated calls).

# Load packages
library(ape)
library(data.table)
library(logging)
library(phangorn)
basicConfig("INFO")

args <- commandArgs(trailingOnly = TRUE)
smpl_assn_file <- args[1]
tree_file <- args[2]

if (!file.exists(smpl_assn_file)) {
    stop("Sample assignment file not found")
}

if (!file.exists(tree_file)) {
    stop("Tree file not found")
}

# Function to print numbers with commas
num = function(x) prettyNum(x, big.mark=",")

# Load variant data
cat("Loading variant tables...\n")
snvs.metadata <-   readRDS("VariantTables_snvs_Metadata.RDS")
snvs.nr <-         readRDS("VariantTables_snvs_NR.RDS")
snvs.nv <-         readRDS("VariantTables_snvs_NV.RDS")
indels.metadata <- readRDS("VariantTables_indels_Metadata.RDS")
indels.nr <-       readRDS("VariantTables_indels_NR.RDS")
indels.nv <-       readRDS("VariantTables_indels_NV.RDS")
for (j in c(1, 3, 4, 5, 7, 8)) {
    snvs.metadata[, j] = as.character(snvs.metadata[, j])
    indels.metadata[, j] = as.character(indels.metadata[, j])
}

# Compute VAF
cat("Computing VAF...\n")
snvs.vaf = snvs.nv / snvs.nr
snvs.vaf[is.na(snvs.vaf)] = 0
indels.vaf = indels.nv / indels.nr
indels.vaf[is.na(indels.vaf)] = 0

# Create tumour and host indices
# TODO: Use external table to define tumour and host samples, not regex
smpl_assn <- fread(smpl_assn_file)
setnames(smpl_assn, c("Sample", "Type"))
if (!(all(smpl_assn$Type %in% c("H", "T")))) {
    logerror("Something is wrong with the sample file. Types other than H|T were supplied in Column 2")
    stop("ERROR - Bad Input")
}
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


# Prepare filter indices
cat("Before filtering:", num(nrow(snvs.metadata)), "SNVs and", num(nrow(indels.metadata)), "indels\n")

# 1) Low-coverage filter
# Select variants with median coverage <10 or >200 across hosts
cat("Running low-coverage filter...\n")
MIN.NR = 10
MAX.NR = 200
median.nr.host.snv = apply(snvs.nr[, hosts], 1, median)
median.nr.host.ind = apply(indels.nr[, hosts], 1, median)
low.nr.snv.idx = (median.nr.host.snv < MIN.NR) | (median.nr.host.snv > MAX.NR)
low.nr.ind.idx = (median.nr.host.ind < MIN.NR) | (median.nr.host.ind > MAX.NR)
cat("Low-coverage filter: discarding", num(sum(low.nr.snv.idx)), "SNVs and",
    num(sum(low.nr.ind.idx)), "indels\n")


# 2) Host filter
# Select variants that occur with 1–4 reads in ≥1 host, and never with ≥5 reads in any sample
cat("Running host filter...\n")
MIN.NV = 1
MAX.NV = 4
low.nv.host.snv.idx = rowSums(snvs.nv[, hosts] >= MIN.NV & snvs.nv[, hosts] <= MAX.NV) > 0 &
    rowSums(snvs.nv > MAX.NV) == 0
low.nv.host.ind.idx = rowSums(indels.nv[, hosts] >= MIN.NV & indels.nv[, hosts] <= MAX.NV) > 0 &
    rowSums(indels.nv > MAX.NV) == 0
cat("Host filter: discarding", num(sum(low.nv.host.snv.idx & !low.nr.snv.idx)), "SNVs and",
    num(sum(low.nv.host.ind.idx & !low.nr.ind.idx)), "indels\n")


# 3) Tumour filter
# Select variants that occur with 1–4 reads in ≥2 tumours, and never with ≥5 reads in any tumour
cat("Running tumour filter...\n")
MIN.NV = 1
MAX.NV = 4
MIN.TUM = 2
low.nv.tum.snv.idx = rowSums(snvs.nv[, tumours] >= MIN.NV & snvs.nv[, tumours] <= MAX.NV) >= MIN.TUM &
    rowSums(snvs.nv[, tumours] > MAX.NV) == 0
low.nv.tum.ind.idx = rowSums(indels.nv[, tumours] >= MIN.NV & indels.nv[, tumours] <= MAX.NV) >= MIN.TUM &
    rowSums(indels.nv[, tumours] > MAX.NV) == 0

# Rescue those SNVs that do not occur in all tumours, and for which the tumours with 1–4 reads
# are phylogenetically related (i.e. the variant represents exactly 1 change in the tree).
cat("Running phylogenetic rescue...\n")

# Construct presence matrix for filter candidates (SNVs ony)
presence = t(snvs.nv[low.nv.tum.snv.idx, tumours] >= MIN.NV)

# Load preliminary tree (by Kevin)
# TODO: replace with pipeline's guide tree
tree = read.tree(tree_file)
stopifnot(identical(sort(tree$tip.label), sort(rownames(presence))))

# Phangorn::pml hangs forever if there are any zero-length branches in the tree,
# so I fix them to an arbitrary small positive value
tree$edge.length[tree$edge.length < 1e-6] <- 1e-6

# (Code below by Kevin Gori: see Scripts/ancestral_state_changes_Kevin.R)
# Compute number of changes per variant.
# Input data format: data frame with rows representing samples and columns representing
# variants. Each value is TRUE if the variant is present, FALSE otherwise.
# Convert presence matrix to phylogenetic data
aln = phangorn::as.phyDat(presence, type="USER", levels=c(F, T))

# Compute tree likelihood
cat(" * Computing tree likelihood...\n")
fit = phangorn::pml(tree, aln)

# ML refine tree [optional] - optimises model parameters
# Set optNni=FALSE to turn off tree topology rearrangements
#phangorn::optim.pml(fit, optNni = TRUE, optBf = TRUE, optQ = TRUE)

# Root the tree - midpoint rooting works in most cases
fit$tree = phangorn::midpoint(fit$tree)
#plot(fit$tree)

# Get ancestral states
#anc.ml = phangorn::ancestral.pml(fit, "ml")        # ML
cat(" * Computing ancestral states...\n")
anc.pars = phangorn::ancestral.pars(fit$tree, aln)  # Parsimony

# Plot
#plotAnc(fit$tree, anc.pars, 1)
#plotAnc(fit$tree, anc.pars, 2)
#plotAnc(fit$tree, anc.pars, 3)  # etc

# Function to traverse the tree from root to tips, adding up the amount of
# state change in each branch. Returns total state change.
count.changes = function(tree, anc.data, site.pattern) {
    total.change = 0
    for (index in rev(ape::postorder(tree))) {
        # An 'edge' consists of a 2-length index of parent and child nodes
        edge = tree$edge[index, ]

        # Look up parent and child ancestral character states using edge index
        parent.state = anc.data[[edge[1]]][site.pattern, ]
        child.state = anc.data[[edge[2]]][site.pattern, ]

        # Compute distance between the parent and child states as a
        # measure of the amount of change on this branch
        # Normalised L1 distance:
        change = sum(abs(parent.state - child.state)) / 2
        # Normalised Euclidean distance:
        #change = sqrt(sum((parent.state - child.state)^2)) / sqrt(2)
        total.change = total.change + change
    }
    total.change
}

# Create data frame mapping site pattern indexes to state change scores
cat(" * Computing state change scores...\n")
scores = data.frame(site.pattern=1:attr(anc.pars, "nr"),
                    pars.score=sapply(1:attr(anc.pars, "nr"),
                                      function(i) count.changes(fit$tree, anc.pars, i)))

# Plot SNV site patterns
cat(" * Plotting SNV site patterns to site_patterns/ folder\n")
dir.create("site_patterns", showWarnings=F)
pdf(paste0("site_patterns/", Sys.Date(), "_TreePatterns_SNVs_1-4reads_1-2changes.pdf"), 6, 9)
for (idx in order(scores$pars.score)) {
    # Plot only for 1-2 changes
    if (scores$pars.score[idx] < 3) {
        plotAnc(tree=fit$tree, data=anc.pars, i=scores$site.pattern[idx], pos=NULL,
                main=sprintf("%.2f changes (parsimony)\n%i variant(s)", scores$pars.score[idx],
                             sum(attr(anc.pars, "index") == scores$site.pattern[idx])))
    }
}
invisible(dev.off())

# Rescue SNVs whose site pattern involves exactly 1 change
phylo.1chg.snv.idx = scores$pars.score[attr(anc.pars, "index")] == 1
low.nv.tum.snv.idx[low.nv.tum.snv.idx][phylo.1chg.snv.idx] = FALSE

cat("Tumour filter: discarding", num(sum(low.nv.tum.snv.idx & !low.nv.host.snv.idx & !low.nr.snv.idx)), "SNVs and",
    num(sum(low.nv.tum.ind.idx & !low.nv.host.ind.idx & !low.nr.ind.idx)), "indels (rescued",
    num(sum(phylo.1chg.snv.idx)), "SNVs)\n")


# 4) Non-chromosomal filter
cat("Running non-chromosomal filter...\n")
non.chrom.snv.idx = !(snvs.metadata$CHROM %in% c(1:38, "X", "Y"))
non.chrom.ind.idx = !(indels.metadata$CHROM %in% c(1:38, "X", "Y"))
cat("Non-chromosomal filter: discarding",
    num(sum(non.chrom.snv.idx & !low.nv.tum.snv.idx & !low.nv.host.snv.idx & !low.nr.snv.idx)), "SNVs and",
    num(sum(non.chrom.ind.idx & !low.nv.tum.ind.idx & !low.nv.host.ind.idx & !low.nr.ind.idx)), "indels\n")


# 5) Duplicate variant filter
cat("Running duplicate variant filter...\n")
dup.snv.idx = duplicated(snvs.metadata[, c("CHROM", "POS", "ALT")])
dup.ind.idx = duplicated(indels.metadata[, c("CHROM", "POS", "ALT")])
cat("Duplicate variant filter: discarding", num(sum(dup.snv.idx)), "duplicated SNVs and",
    num(sum(dup.ind.idx)), "duplicated indels\n")


# Apply all filters
cat("Applying filters...\n")
filter.snv.idx = low.nr.snv.idx | low.nv.host.snv.idx | low.nv.tum.snv.idx | non.chrom.snv.idx | dup.snv.idx
filter.ind.idx = low.nr.ind.idx | low.nv.host.ind.idx | low.nv.tum.ind.idx | non.chrom.ind.idx | dup.ind.idx

snvs.metadata = snvs.metadata[!filter.snv.idx, ]
snvs.nr = snvs.nr[!filter.snv.idx, ]
snvs.nv = snvs.nv[!filter.snv.idx, ]
snvs.vaf = snvs.vaf[!filter.snv.idx, ]
indels.metadata = indels.metadata[!filter.ind.idx, ]
indels.nr = indels.nr[!filter.ind.idx, ]
indels.nv = indels.nv[!filter.ind.idx, ]
indels.vaf = indels.vaf[!filter.ind.idx, ]
cat("After filtering:", num(nrow(snvs.metadata)), "SNVs and", num(nrow(indels.metadata)), "indels\n")


# Save objects
cat("Saving tables to data/ directory...\n")
dir.create("data")
save(snvs.metadata,   file = "data/VariantTables_Filtered_snvs.metadata.RData")
save(snvs.nr,         file = "data/VariantTables_Filtered_snvs.nr.RData")
save(snvs.nv,         file = "data/VariantTables_Filtered_snvs.nv.RData")
save(snvs.vaf,        file = "data/VariantTables_Filtered_snvs.vaf.RData")
save(indels.metadata, file = "data/VariantTables_Filtered_indels.metadata.RData")
save(indels.nr,       file = "data/VariantTables_Filtered_indels.nr.RData")
save(indels.nv,       file = "data/VariantTables_Filtered_indels.nv.RData")
save(indels.vaf,      file = "data/VariantTables_Filtered_indels.vaf.RData")
