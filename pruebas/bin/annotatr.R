#!/usr/bin/env Rscript

# LOAD REQUIRED LIBRARIES
library(GenomicRanges)
library(annotatr)
library(WriteXLS)

# Get command line arguments
args <- commandArgs(trailingOnly = TRUE)
sample <- args[1]
allsites_file <- args[2] 
ref_genome <- args[3]  # hg38, hg19, mm10, etc.

# Función para anotar con annotatr - ACEPTA NCBI
annotate_with_annotatr <- function(granges_obj, genome_name) {
  message("Anotando con annotatr para genoma: ", genome_name)

  # Convertir nombre genoma a formato annotatr
  annotatr_genome <- switch(genome_name,
    "hg38" = "hg38",
    "hg19" = "hg19", 
    "mm10" = "mm10",
    "mm9" = "mm9",
    "rn6" = "rn6",
    stop("Genoma no soportado por annotatr: ", genome_name)
  