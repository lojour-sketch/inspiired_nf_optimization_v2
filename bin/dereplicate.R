#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)

sample <- args[1]
read1 <- args[2] 
read2 <- args[3]

library(ShortRead)

cat("[" , Sys.time(), "] Reading R1...\n")
r1 <- readFastq(read1)
cat("[" , Sys.time(), "] R1 loaded:", length(r1), "reads\n")

cat("[" , Sys.time(), "] Reading R2...\n")
r2 <- readFastq(read2)
cat("[" , Sys.time(), "] R2 loaded:", length(r2), "reads\n")

cat("[" , Sys.time(), "] Finding unique sequences...\n")
# Get sequences as character for matching
r1_seqs <- as.character(sread(r1))
r2_seqs <- as.character(sread(r2))

# Debug: Check for empty sequences
cat("DEBUG: Checking for empty sequences in R1...\n")
empty_seqs <- which(width(sread(r1)) == 0)
if(length(empty_seqs) > 0) {
    cat("WARNING: Found", length(empty_seqs), "empty sequences in R1\n")
    print(empty_seqs[1:3])
}
empty_seqs_r2 <- which(width(sread(r2)) == 0)
if(length(empty_seqs_r2) > 0) {
    cat("WARNING: Found", length(empty_seqs_r2), "empty sequences in R2\n")
    print(empty_seqs_r2[1:3])
}
# Remove empty sequences before proceeding
if(length(empty_seqs) > 0) {
    r1 <- r1[-empty_seqs]
    cat("Removed", length(empty_seqs), "empty sequences from R1\n")
}
if(length(empty_seqs_r2) > 0) {
    r2 <- r2[-empty_seqs_r2]
    cat("Removed", length(empty_seqs_r2), "empty sequences from R2\n")
}



# Find first occurrence of each unique sequence
r1_unique_idx <- which(!duplicated(r1_seqs))
r2_unique_idx <- which(!duplicated(r2_seqs))

# Keep the original objects with original IDs and qualities
r1u <- r1[r1_unique_idx]
r2u <- r2[r2_unique_idx]

cat("[" , Sys.time(), "] R1 unique:", length(r1u), "sequences\n")
cat("[" , Sys.time(), "] R2 unique:", length(r2u), "sequences\n")

cat("[" , Sys.time(), "] Writing unique files...\n")
writeFastq(r1u, paste0(sample, ".R1_unique_by_sequence.fastq"), compress = FALSE)
writeFastq(r2u, paste0(sample, ".R2_unique_by_sequence.fastq"), compress = FALSE)

cat("[" , Sys.time(), "] Creating keys dataframe...\n")
keys <- data.frame(name = as.character(id(r1)),
                   R1 = match(r1_seqs, r1_seqs[r1_unique_idx]),
                   R2 = match(r2_seqs, r2_seqs[r2_unique_idx]))
cat("[" , Sys.time(), "] Keys dataframe created\n")

cat("[" , Sys.time(), "] Saving RData...\n")
save(keys, file = paste0(sample, "_keys.RData"))
cat("[" , Sys.time(), "] Done!\n")