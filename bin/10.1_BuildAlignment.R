#!/usr/bin/env Rscript
# Analysis of CTVT whole-genome variants
# Step 10.1. Build alignment of tumour-only SNVs for tree inference
# Adrian Baez-Ortega, 2018

# As an input to RAxML, we generate an alignment of the reference sequence and all
# the tumours, consisting of the concatenation of all somatic SNVs in positions
# with sufficient read coverage. The SNVs will be separated into coding and
# noncoding, such that two partitions can be defined. Allele imputation process
# is similar to that previously used for the CTVT exomes:
# We calculate the power to call a variant if it is present at each position. If
# we have insufficient power, we impute the IUPAC character representing ambiguity
# between the REF and ALT alleles. Only if there is sufficient power we impute the
# ALT or REF allele. We define sufficient power as power to sequence 10 supporting
# reads. For each position in each sample, we estimate the number of reads we would
# expect in a heterozygous situation using coverage and purity. If <10, then we
# impute the IUPAC ambiguity character. If ≥10, we impute the ALT allele for that
# variant if it is supported by ≥5 reads; otherwise, we impute it as REF.


# Minimum read thresholds for non-ambiguity and support
MIN.EXPECTED = 10
MIN.SUPPORT = 5

# IUPAC ambiguity character matrix
IUPAC = matrix(c("-", "M", "R", "W",
                 "M", "-", "S", "Y",
                 "R", "S", "-", "K",
                 "W", "Y", "K", "-"),
               nrow=4, byrow=T,
               dimnames=list(c("A", "C", "G", "T"), c("A", "C", "G", "T")))

# Function to print numbers with commas
num = function(x) prettyNum(x, big.mark=",")

# Load variant data
load("TumourPurity.RData")
load("VariantTables_Indices_Filt.RData")
snvs.nr.tonly    <- readRDS("VariantTables_Split_Filt_snvs_NR.tonly.RDS")
snvs.nv.tonly    <- readRDS("VariantTables_Split_Filt_snvs_NV.tonly.RDS")
annot.snvs.tonly <- readRDS("AnnotationTables_annot.snvs.tonly.RDS")


# Define variant order: coding (CDS) SNVs first, non-coding SNVs after
coding.idx = annot.snvs.tonly$Codons != "-"
stopifnot(!any(is.na(coding.idx)))
order.idx = c(which(coding.idx), which(!coding.idx))

# Keep track of last coding SNV in alignment
last.coding.idx = sum(coding.idx)


# Initialise alignment
cat("Imputing", num(sum(coding.idx)), "coding and", num(sum(!coding.idx)),
    "non-coding SNVs on", sum(tumours), "tumours...\n")

imputed.alleles = matrix(rep(annot.snvs.tonly$Ref[order.idx],
                             each=sum(tumours) + 1),
                         nrow=sum(tumours) + 1,
                         dimnames=list(c("REF", samples[tumours]), NULL))

alt.alleles = matrix(rep(annot.snvs.tonly$Alt[order.idx],
                         each=sum(tumours) + 1),
                     nrow=sum(tumours) + 1,
                     dimnames=list(c("REF", samples[tumours]), NULL))

metaidx = annot.snvs.tonly$Meta_index[order.idx]


# Check that matrix*vector product works as expected
stopifnot(identical(matrix(1:12, 3) * (1:3), 
                    matrix(as.integer(c(1,4,9,4,10,18,7,16,27,10,22,36)), 3)))

# Calculate expected number of reads for a HET variant
expected.reads = round(t(snvs.nr.tonly[metaidx, tumours]) * purity / 2)

# Identify variants with insufficient expected reads in each tumour
min.expected.idx = expected.reads >= MIN.EXPECTED

# Identify variants with sufficient read support
min.support.idx = t(snvs.nv.tonly[metaidx, tumours]) >= MIN.SUPPORT


# Impute IUPAC character in tumours with insufficient expected reads
insuf.exp.ref = imputed.alleles[-1, ][!min.expected.idx]
insuf.exp.alt = alt.alleles[-1, ][!min.expected.idx]
imputed.alleles[-1, ][!min.expected.idx] = sapply(1:length(insuf.exp.ref), function(i) {
    IUPAC[insuf.exp.ref[i], insuf.exp.alt[i]]
})

# Impute ALT allele in tumours with sufficient variant support
imputed.alleles[-1, ][min.expected.idx & min.support.idx] = 
    alt.alleles[-1, ][min.expected.idx & min.support.idx]


# Discard variants for which no allele is ALT
discard.idx = colSums(min.expected.idx & min.support.idx) == 0

cat("Discarding", num(sum(discard.idx)), "uninformative sites (SNVs)...\n")
imputed.alleles = imputed.alleles[, !discard.idx]
last.coding.idx = last.coding.idx - sum(which(discard.idx) <= last.coding.idx)

cat("Final alignment has ", nrow(imputed.alleles), " samples and ", num(ncol(imputed.alleles)), 
    " sites (", num(last.coding.idx), " coding sites)\n\n", sep="")


# Output alignment in PHYLIP format
cat("Writing alignment (PHYLIP) and last coding site position to data/ directory...\n")
dir.create("data")
write(dim(imputed.alleles), file="data/TumourAlignment.phylip")
for (i in 1:nrow(imputed.alleles)) {
    write(c(rownames(imputed.alleles)[i], paste(imputed.alleles[i, ], collapse="")),
          file="data/TumourAlignment.phylip", ncolumns=2, append=T)
}

# Save index of last coding position
write(last.coding.idx, file="data/TumourAlignment_LastCodingSite.txt")

cat(sprintf("ASC_DNAX, coding = 1 - %d\n", last.coding.idx),
    sprintf("ASC_DNAX, noncoding = %d - %d\n", last.coding.idx + 1, ncol(imputed.alleles)),
    file="data/partitions.txt", sep="", append=T)
