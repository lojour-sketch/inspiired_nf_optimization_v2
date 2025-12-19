#!/usr/bin/env Rscript

#load arguments
args <- commandArgs(trailingOnly = TRUE)
fq1 <- args[1] #tsv
fq2 <- args[2] #tsv
sample <- args[3] #character

library(ShortRead)
library(FastqCleaner)

# Read FASTQ
message("Reading FASTQ files...")
fq1 <- readFastq(fq1)
fq2 <- readFastq(fq2)

#we homogeneize the names of the reads, so that they are the same in both R1 and R2 files
new_ids <- sapply(
  as.character(id(fq2)),        # convertir BStringSet a character
  function(z) paste0(sample, "%", strsplit(z, "-")[[1]][2])
)

message("Homogeneizing read names in R1 and R2 files...")
fq1 <- ShortReadQ(sread(fq1), quality(fq1), BStringSet(new_ids))
fq2 <- ShortReadQ(sread(fq2), quality(fq2), BStringSet(new_ids))


# Extract read IDs as character vectors
ids1 <- as.character(id(fq1))
ids2 <- as.character(id(fq2))

#we do intersection, as in inspiired. 
message("Intersecting R1 and R2 files...")
common <- intersect(ids1, ids2)
fq1 <- fq1[ ids1 %in% common ]
fq2 <- fq2[ ids2 %in% common ]

# Get unique sequences from filtered reads, to do the matching in the keys dataframe
message("Extracting unique sequences from filtered reads...")
fq1.u <- unique(sread(fq1))
fq2.u <- unique(sread(fq2))

message("Assigning numbers to unique sequences...")
names(fq1.u) <- as.character(seq_along(fq1.u))
names(fq2.u) <- as.character(seq_along(fq2.u))

# Build keys with match, as in original logic
message("Building keys dataframe...")
keys <- data.frame(
  R2 = match(as.character(sread(fq2)), as.character(fq2.u)),
  R1 = match(as.character(sread(fq1)), as.character(fq1.u)),
  names = as.character(id(fq2))
)
keys$readPairKey <- paste0(keys$R1, "_", keys$R2)

save(keys, file = paste0(sample, "_keys.RData"))

#now we will create the fastq files with numerical headers, to be used in the postprocessing with the keys data. 
message("deduplicating and writing FASTQ files with numerical headers...")
# Deduplicate sequences to match keys
fq1_dedup <- unique_filter(fq1)
fq2_dedup <- unique_filter(fq2)

# message("Intersecting R1 and R2 files...")
# common <- intersect(ids1, ids2)
# fq1 <- fq1[ ids1 %in% common ]
# fq2 <- fq2[ ids2 %in% common ]

# # Sanity check: lengths must match after deduplication
# stopifnot(length(fq1_dedup) == length(fq2_dedup))

# Assign numeric IDs (1..n) to deduplicated reads
num_ids1 <- as.character(seq_along(fq1_dedup))
num_ids2 <- as.character(seq_along(fq2_dedup))

message("Changing to numeric read names in R1 and R2 files...")
fq1_dedup <- ShortReadQ(sread(fq1_dedup), quality(fq1_dedup), BStringSet(as.character(num_ids1)))
fq2_dedup <- ShortReadQ(sread(fq2_dedup), quality(fq2_dedup), BStringSet(as.character(num_ids2)))

message("Printing first 4 lines of R1 and R2 files...")
print(head(id(fq1_dedup), 4))
print(head(id(fq2_dedup), 4))

message("Writing numeric FASTQ files...")
writeFastq(fq1_dedup, paste0(sample, "_R1_numeric.fastq.gz"), compressed = TRUE)
writeFastq(fq2_dedup, paste0(sample, "_R2_numeric.fastq.gz"), compressed = TRUE)

message("Numeric FASTQ files ready for downstream postprocessing.")

# Optional: create a mapping from numeric IDs back to original reads
# Useful if you need to expand results later
# mapping <- data.frame(
#   originalID_R1 = as.character(id(fq1_dedup)),
#   originalID_R2 = as.character(id(fq2_dedup))
# )
# save(mapping, file = paste0(sample, "_numericID_mapping.RData"))

# message("Mapping object saved for potential re-expansion of duplicates. However, INSPIIRED doesn't seem to do that even if it says it does.")