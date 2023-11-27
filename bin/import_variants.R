#!/usr/bin/env Rscript

library(logging)
basicConfig("INFO")

filename = commandArgs(trailing = TRUE)[1]
outfile = commandArgs(trailing = TRUE)[2]

if (!file.exists(filename)) {
    logerror(paste("File", filename, "not found"))
    stop("ERROR - file error")
}

check_names <-  grepl("Metadata", filename)
tbl = read.table(filename, header = TRUE, check.names = check_names)
save(tbl, file = outfile)
rm(tbl)

