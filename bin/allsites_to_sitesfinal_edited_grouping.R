#!/usr/bin/env Rscript

# ------------------------------------------------------------------
# Author: Libe Renteria Aizpurua
# Date: 2026-01-07 
#
# This script groups all insertion sites of allSites that have identical coordinates and strand.
# It also counts the number of times each group appears in the allSites file.
# It also creates a list of indexes of the reads where that unique siteid appears.
# This way we can see the level of clonal expansion of the samples.
#
# ------------------------------------------------------------------


#load libraries
library(GenomicRanges)
library(hiReadsProcessor)
library(parallel)
library(rtracklayer)
library(readr)
library(dplyr)
library(S4Vectors) 


#import command line arguments
args <- commandArgs(trailingOnly = TRUE)
sample <- args[1]
allsites_file <- args[2]

dereplicateSites_edited <- function(uniqueReads){
  if (length(uniqueReads) == 0) return("No sites to dereplicate")
  print(paste0("dereplicating ", length(uniqueReads), " sites"))
  
  #obtain an unique identificator per start, end and strand
  siteIDs <- paste0(seqnames(uniqueReads), ":", start(uniqueReads), "-", end(uniqueReads), ":", strand(uniqueReads))

  #count occurrences of each siteID
  siteIDs.counts <- table(siteIDs)

  #get only one of each siteID
  siteIDs.unique <- names(siteIDs.counts)

  # revmap: list of indexes of the reads where that unique siteid appears
  #with which we obtain the index of the first match of each element in siteIDs.unique to siteIDs, this is, the index of the reads where that unique siteid appears
  revmap_idx <- lapply(siteIDs.unique, function(id) which(siteIDs == id))
  revmap_idx <- IRanges::IntegerList(revmap_idx) 

  #the match function returns the index of the first match of each element in siteIDs.unique to siteIDs, so we only keep the ranges in those indices
  dereplicatedSites <- uniqueReads[match(siteIDs.unique, siteIDs)]

  #add counts and revmap as metadata. The rest of the metadata is already in the GRanges object
  mcols(dereplicatedSites)$revmap <- revmap_idx
  mcols(dereplicatedSites)$counts <- as.integer(siteIDs.counts)
  
  dereplicatedSites
}

#first we convert our allSites.tsv file to a dataframe
colClasses <- c("character", "integer", "integer", "character", "list", "integer", "character", "character")
df <- read.table(allsites_file, header=TRUE, sep="\t", stringsAsFactors=FALSE, row.names=8, colClasses=colClasses)

allSites_gr <- makeGRangesFromDataFrame(df, keep.extra.columns=TRUE, seqnames.field="seqnames", start.field="start", end.field="end", strand.field="strand")

sites.final <- dereplicateSites_edited(allSites_gr)

saveRDS(allSites_gr, file=paste0(sample, "_allsites.rds"))
saveRDS(sites.final, file=paste0(sample, "_sitesfinal.rds"))
