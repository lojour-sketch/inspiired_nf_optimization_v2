#!/usr/bin/env Rscript

# Script para estimar uso real de memoria

args <- commandArgs(trailingOnly = TRUE)
if(length(args) < 1) {
    cat("Usage: estimate_memory_usage.R <fastq.gz file>\n")
    quit(status=1)
}

read1 <- args[1]

library(ShortRead)
library(Biostrings)
library(pwalign)

cat("=== MEMORY USAGE ESTIMATION ===\n\n")

# Tamaño del archivo
file_size_gb <- file.size(read1) / (1024^3)
cat(sprintf("Input file size: %.2f GB (compressed)\n", file_size_gb))

# Estimar tamaño descomprimido (FASTQ comprime ~3-4x)
estimated_uncompressed <- file_size_gb * 3.5
cat(sprintf("Estimated uncompressed: %.2f GB\n", estimated_uncompressed))

cat("\nReading file and measuring actual memory usage...\n")

# Medir memoria antes
gc()
mem_before <- sum(gc()[,2])

# Leer archivo
fq1 <- readFastq(read1)
reads1 <- sread(fq1)
qual1 <- quality(fq1)

# Medir memoria después
gc()
mem_after <- sum(gc()[,2])
mem_used_mb <- mem_after - mem_before

cat(sprintf("\n=== RESULTS ===\n"))
cat(sprintf("Number of reads: %d\n", length(fq1)))
cat(sprintf("Memory used for loading: %.2f MB (%.2f GB)\n", 
            mem_used_mb, mem_used_mb/1024))

# Test alignment memory
cat("\nTesting alignment memory with 1000 reads...\n")
marker <- "TGTGGAAAGGACGAAACACCG"
marker_rc <- reverseComplement(DNAString(marker))
submat <- pwalign::nucleotideSubstitutionMatrix(match=1, mismatch=0, baseOnly=TRUE)

gc()
mem_before_align <- sum(gc()[,2])

test_reads <- reads1[1:1000]
tmp <- pwalign::pairwiseAlignment(
    pattern=test_reads,
    subject=marker_rc,
    substitutionMatrix=submat,
    gapOpening=0,
    gapExtension=1,
    type="overlap"
)

gc()
mem_after_align <- sum(gc()[,2])
align_mem_per_1k <- (mem_after_align - mem_before_align)
align_mem_total <- (align_mem_per_1k / 1000) * length(reads1)

cat(sprintf("Alignment memory (1k reads): %.2f MB\n", align_mem_per_1k))
cat(sprintf("Estimated alignment memory (all): %.2f MB (%.2f GB)\n", 
            align_mem_total, align_mem_total/1024))

# Peak memory estimate
peak_memory_gb <- (mem_used_mb + align_mem_total) / 1024 * 2.5  # 2.5x safety factor

cat("\n=== RECOMMENDATION ===\n")
cat(sprintf("Estimated peak memory usage: %.2f GB\n", peak_memory_gb))
cat(sprintf("Recommended memory setting: %.0f GB\n", ceiling(peak_memory_gb * 1.3)))
cat("\nNote: This includes 30%% safety margin\n")

if(peak_memory_gb < 50) {
    cat("You're requesting 100GB but probably only need ~", 
        ceiling(peak_memory_gb * 1.3), "GB\n", sep="")
    cat("Consider reducing memory request to improve parallelization!\n")
}