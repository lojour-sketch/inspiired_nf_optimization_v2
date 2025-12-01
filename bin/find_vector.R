#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)

# Parse arguments
meta <- args[1]
read1 <- args[2] 
read2 <- args[3]
primerltr <- args[4]
vector <- args[5]
globalIdentity <- 0.75

library(Biostrings)
library(ShortRead)

## we will run the alignment with minimap. However, minimap calculates the identity percentage differently from BLAT, because minimap takes into account gaps to calculate the query length
## so we will calculate the identity the BLAT way (as inspiired) calculating matches / qSize from minimaps output
## then we will add this new identity to the output PAF file as the 13th column
calculateidentity <- function(query_file, subject_file, min_identity) {
    # Run minimap2
    paf_file <- tempfile(pattern = "minimap", fileext = ".paf")
    system2("minimap2", 
            args = c("-x", "sr", "-N", "100", "-p", "0.1", "--eqx",
                    subject_file, query_file),
            stdout = paf_file, stderr = FALSE)
        
    if(file.exists(paf_file) && file.info(paf_file)$size > 0) {
        paf <- read.table(paf_file, sep="\t", fill=TRUE, stringsAsFactors=FALSE)
        # we create a data frame that saves all the info
        result <- data.frame (
            qName = paf[,1],
            qSize = paf[,2], 
            qStart = paf[,3],
            matches = paf[,10],
            tStart = paf[,8],
            identity = paf[,10] / paf[,2],  # BLAT-style identity
            stringsAsFactors = FALSE
            )
        file.remove(paf_file)
        return(result)
    } else {
        return(data.frame())
    }
}

## we align the reads
hits.v.1 <- try(calculateidentity(read1, vector, globalIdentity)) 
hits.v.2 <- try(calculateidentity(read2, vector, globalIdentity))

## we filter exactly as inspiired - using base R instead of dplyr
if (nrow(hits.v.1) > 0) {
    hits.v.1 <- hits.v.1[hits.v.1$matches > globalIdentity * hits.v.1$qSize & hits.v.1$qStart <= 5, ]
}
if (nrow(hits.v.2) > 0) {
    hits.v.2 <- hits.v.2[hits.v.2$matches > globalIdentity * hits.v.2$qSize, ]
}

## Extract common base name from headers (remove the " 1:N:0" / " 2:N:0" part)
get_base_name <- function(header) {
    sapply(strsplit(header, " "), function(x) x[1])
}

## In the merge - use base names for matching
if (nrow(hits.v.1) > 0) {
    hits.v.1$baseName <- get_base_name(hits.v.1$qName)
}
if (nrow(hits.v.2) > 0) {
    hits.v.2$baseName <- get_base_name(hits.v.2$qName)
}

## now we merge the dataframes to get the reads that have vector in both pairs. we will only remove these    
hits.v <- try({
    if (nrow(hits.v.1) > 0 && nrow(hits.v.2) > 0) {
        merge(hits.v.1[, c("baseName", "tStart")],
              hits.v.2[, c("baseName", "tStart")],
              by="baseName")
    } else {
        data.frame(baseName = character(0))
    }
}, silent = TRUE)

## as multihits can happen, we take only unique names
if (inherits(hits.v, "try-error") || nrow(hits.v) == 0) {
    vqName_base <- character(0)
} else {
    vqName_base <- unique(hits.v$baseName) 
}
    
#now we remove reads from both files. first we read them and get their basenames
reads1 <- readFastq(read1, with.qualities=TRUE, format="fastq")
reads2 <- readFastq(read2, with.qualities=TRUE, format="fastq")

# Extract headers safely
get_base_name_safe <- function(shortread_obj) {
    headers <- as.character(id(shortread_obj))
    # Extract base names (part before first space)
    sapply(strsplit(headers, " "), function(x) x[1])
}

reads.1_base <- get_base_name_safe(reads1)
reads.2_base <- get_base_name_safe(reads2)

## Remove vector reads from both files
## vqName_base contains the base names of the headers that belong to the pairs that both have vectors 
reads.1_clean <- reads1[!reads.1_base %in% vqName_base]
reads.2_clean <- reads2[!reads.2_base %in% vqName_base]

output_file1 <- paste0(meta, ".R1_vector_removed.fastq.gz")
output_file2 <- paste0(meta, ".R2_vector_removed.fastq.gz")

if (file.exists(output_file1)) file.remove(output_file1)
if (file.exists(output_file2)) file.remove(output_file2)

# Write FASTQ files properly using ShortRead::writeFastq
if (length(reads.1_clean) > 0) {
    writeFastq(reads.1_clean, output_file1, compress = TRUE)
} else {
    # Create empty file if no reads
    file.create(output_file1)
}

if (length(reads.2_clean) > 0) {
    writeFastq(reads.2_clean, output_file2, compress = TRUE)
} else {
    # Create empty file if no reads
    file.create(output_file2)
}

cat("Successfully processed", length(reads.1_clean), "R1 reads and", length(reads.2_clean), "R2 reads\n")
cat("Removed", length(reads1) - length(reads.1_clean), "R1 reads and", length(reads2) - length(reads.2_clean), "R2 reads containing vector\n")