#!/usr/bin/env Rscript
#instead of doing dereplication (to see repres. starts and order), standardization and posterior dereplication...
#we can do an initial dereplication and standardize the groups, so we don't have to standardize each read and reorder them
#we prefer this script because in the original one the object after standardization was not matching the original length

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

dereplicateSites_Libe <- function(uniqueReads){
  print(paste0("dereplicating ", length(uniqueReads), " sites"))
  #do the dereplication, but loose the coordinates
  sites.reduced <- reduce(flank(uniqueReads, -5, both=TRUE), with.revmap=T, ignore.strand=FALSE)
  sites.reduced$counts <- sapply(sites.reduced$revmap, length)
  
  #order the unique sites as described by revmap
  dereplicatedSites <- uniqueReads[unlist(sites.reduced$revmap)]
  
  #if no sites are present, skip this step - keep doing the rest to provide a
  #similar output to a successful dereplication
  if(length(uniqueReads)>0){
    #split the unique sites as described by revmap (sites.reduced$counts came from revmap above)
    dereplicatedSites <- split(dereplicatedSites, Rle(values=seq(length(sites.reduced)), lengths=sites.reduced$counts))
  }
  
  #do the standardization - this will pick a single starting position and
  #choose the longest fragment as ending position
  dereplicatedSites <- unlist(reduce(flank(dereplicatedSites, -5, both=TRUE)))
  mcols(dereplicatedSites) <- mcols(sites.reduced)
  
  dereplicatedSites
}

#first we load our allSites.RData
object <- load(allsites_file)
allSites_gr <- get(object)

sites.final <- dereplicateSites_Libe(allSites_gr)

saveRDS(allSites_gr, file=paste0(sample, "_allsites_nostandard.rds"))
saveRDS(sites.final, file=paste0(sample, "_sitesfinal.rds"))
