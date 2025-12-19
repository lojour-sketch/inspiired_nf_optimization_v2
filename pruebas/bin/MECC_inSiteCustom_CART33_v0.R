#!/usr/bin/env Rscript
#===============================================================================
#' Author: Maria E. Calleja
#' Date: 2023/MAYJUN
#===============================================================================
## PACKAGES.
library(RMySQL, quietly=TRUE, verbose=FALSE)
library(dplyr, quietly=TRUE, verbose=FALSE)
library(DBI, quietly=TRUE, verbose=FALSE)
library(yaml, quietly=TRUE, verbose=FALSE)
#library()
library(clusterProfiler)
library(org.Hs.eg.db)
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(org.Hs.eg.db)
library(BSgenome.Hsapiens.UCSC.hg19)
library(ggupset)
library(ggimage)

options(stringsAsFactors=F, useFancyQuotes=F)

#===============================================================================
## GLOBAL VARIABLES.
PARENT_DIR<-("/home/mecc/clustermecc/INSPIIRED/CD33_worked") 
"/home/mecc/clustermecc/"
#"INSPIIRED_2022/SHARED/HL_Samples_fangorn/HL_Samples/"
PROJECT_HD68_DIR<-file.path(PARENT_DIR,"CD33_CD68_march")
PROJECT_HD70_DIR<-file.path(PARENT_DIR,"CD33_CD70_march")
HD68_FASTQ_DIR <- file.path(PROJECT_HD68_DIR,"Data/demultiplexedReps")
HD70_FASTQ_DIR <- file.path(PROJECT_HD70_DIR,"Data/demultiplexedReps")

EDITED_DIR<-file.path(PARENT_DIR,"src")
FIGURES_DIR <-file.path(PARENT_DIR,"figures")
RESULTS_DIR <-file.path(PARENT_DIR,"results")

library(BiocParallel)
CORES <- 6
register(SnowParam(CORES))
#===============================================================================
### FROM: 
### "/home/mecc/Apps/INSPIIRED/components/intSiteUploader/intSiteUploader.R"

workingDir <- PARENT_DIR
# Determine analysis directory and change into it.
if( is.na(workingDir) ) workingDir <- "."
workingDir <- normalizePath(workingDir, mustWork=TRUE)
setwd(workingDir)
list.files()
message("Changed to directory: ", workingDir)

dirs <- dir(".", pattern="^CD33")

fileName <- file.path(dirs, "completeMetadata.RData")
testFileName <- sapply(seq_along(dirs), function(x) gsub("_march", "" , dirs[x]))
#https://stackoverflow.com/questions/18751361/how-to-loop-over-files-in-different-directories

for (i in seq(along=fileName)){
  # Get the sample information.
  stopifnot(file.exists(fileName[i]))
  metadata <- get(load(fileName[i]))
  metadata <- subset(metadata, select=c("alias", "gender", "refGenome"))
  names(metadata) <- c("sampleName", "gender", "refGenome")
  print(head(metadata))
  
  message("\nProcessing: ",fileName[i])
  #assign(paste0("Sites_for_",file), sites ,envir = .GlobalEnv)
  
  currentMaxSiteID<-100000000
  currentMaxMultihitID<-100000000
  
  # for(j in seq(nrow(metadata))){
  #   file <- metadata[j,"sampleName"]
  #   message("\nProcessing: ", file)}
  # 
  # for(j in seq(nrow(metadata))){
  #   file <- metadata[j,"sampleName"]
  #   message("\nProcessing: ", file)
  #   if(all(file.exists(file.path(dirs[i],file, "sites.final.RData"), file.path(dirs[i],file, "allSites.RData")))){
  #     load(file.path(dirs[i],file, "sites.final.RData"))
  #     load(file.path(dirs[i],file, "allSites.RData"))} }
  directorio <- dirs[i]
  
  for(j in seq(nrow(metadata))){
    file <- metadata[j,"sampleName"]
    message("\nProcessing: ", file)
    if(all(file.exists(file.path(directorio,file, "sites.final.RData"), file.path(directorio,file, "allSites.RData")))){
      load(file.path(directorio,file, "sites.final.RData"))
      load(file.path(directorio,file, "allSites.RData"))
      if(length(sites.final)>0){
        ##sites.final won't exist if there aren't sites, thus no need to check if sites.final has sites in it
        sites <- data.frame(
          "siteID"=seq(length(sites.final))+currentMaxSiteID,
          "sampleID"=metadata[j,"sampleName"],
          "position"=start(flank(sites.final, -1, start=T)),
          "chr"=as.character(seqnames(sites.final)),
          "strand"=as.character(strand(sites.final)) )
        
        ## change to the right class as database
        sites$siteID <- as(sites$siteID, "integer")
        sites$sampleID <- as(sites$sampleID, "character")
        sites$position <- as(sites$position, "integer")
        sites$chr <- as(sites$chr, "character")
        sites$strand <- as(sites$strand, "character")
        
        ##Newer versions of intSiteCaller return allSites in the order dictated by
        ##sites.final.  This line allows import of 'legacy' output
        allSites <- allSites[unlist(sites.final$revmap)]
        
        ##could do the next three statements with aggregate, but this method is emperically 2x faster
        pcrBreakpoints <- sort(paste0(as.integer(Rle(sites$siteID, sapply(sites.final$revmap, length))),
                                      ".",
                                      start(flank(allSites, -1, start=F))))
        
        condensedPCRBreakpoints <- strsplit(unique(pcrBreakpoints), "\\.")
        
        pcrBreakpoints <- data.frame("siteID"=sapply(condensedPCRBreakpoints, "[[", 1),
                                     "breakpoint"=sapply(condensedPCRBreakpoints, "[[", 2),
                                     "count"=runLength(Rle(match(pcrBreakpoints, unique(pcrBreakpoints)))))
        
        # change to the right class as database
        pcrBreakpoints$siteID <- as(pcrBreakpoints$siteID, "integer")
        pcrBreakpoints$breakpoint <- as(pcrBreakpoints$breakpoint, "integer")
        pcrBreakpoints$count <- as(pcrBreakpoints$count, "integer")
        
        # MINE: instead of loading it into the database, Ill read them here. 
        assign(paste0("Sites_for_",file), sites ,envir = .GlobalEnv)
        assign(paste0("PCRbreakpoints_for_",file), pcrBreakpoints ,envir = .GlobalEnv)
        
        # load table sites            
        message("Loading sites: ", nrow(sites), " entries")
        ### stopifnot( dbWriteTable(dbConn, "sites", sites, append=T, row.names=F) )
        #null <- apply(sites, 1, insertSQL, dbConn=dbConn, colNames=colnames(sites), table='sites')            
        
        # load table pcrbreakpoints
        message("Loading pcrbreakpoints: ", nrow(pcrBreakpoints), " entries")
        ### stopifnot( dbWriteTable(dbConn, "pcrbreakpoints", pcrBreakpoints, append=T, row.names=F) )
        #null <- apply(pcrBreakpoints, 1, insertSQL, dbConn=dbConn, colNames=colnames(pcrBreakpoints), table='pcrbreakpoints')            
        
        newMaxSiteID = currentMaxSiteID + nrow(sites)
        currentMaxSiteID <- newMaxSiteID
        
        peak <- sites.final
        cromosomas = c(paste0("chr",1:23), "chrX")
        txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
        peakAnno <- annotatePeak(peak, tssRegion=c(-3000, 3000),
                                 TxDb=txdb, annoDb="org.Hs.eg.db")
        assign(paste0("peakAnno_",file), peakAnno ,envir = .GlobalEnv)
        peakAnno.dfr <- as.data.frame(as.data.frame(peakAnno))
        
        #WriteXLS::WriteXLS(peakAnno.dfr, paste0(FIGURES_DIR,"/annotated_", testFileName[i], file,'.xlsx'))
        WriteXLS::WriteXLS(peakAnno.dfr, paste0(FIGURES_DIR,"/annotated_", file,'.xlsx'))
        #write.table(peakAnno.dfr, paste0(FIGURES_DIR,"/annotated_", file,'.txt'), sep="\t", quote=F, col.names=T, row.names=F)
        #write.csv(peakAnno.dfr, paste0(FIGURES_DIR,"/annotated_", file,'.csv'))
        #write.csv(peakAnno.dfr, paste0(FIGURES_DIR,"/annotated_", file,'.txt'))
        
        pdf(paste0(FIGURES_DIR,"/fig_annot_", file,'.pdf'))
        print(covplot(peak, weightCol = "counts", chrs= cromosomas, title = paste0("Insertion Sites over Chromosomes sample: ", file)))
        print(plotAnnoPie(peakAnno))
        print(plotAnnoBar(peakAnno))
        print(upsetplot(peakAnno))
        print(vennpie(peakAnno))
        print(upsetplot(peakAnno, vennpie=TRUE))
        dev.off()
      }
    }
    
    if(file.exists(file.path(directorio,file, "multihitData.RData"))){
      load(file.path(directorio,file, "multihitData.RData"))
      
      if(length(multihitData[[1]])>0){
        multihitPositions <- multihitData[[2]]
        multihitLengths <- multihitData[[3]]
        stopifnot(length(multihitPositions)==length(multihitLengths))
        multihitPositions <- data.frame(
          "multihitID"=rep(seq(length(multihitPositions))+currentMaxMultihitID,
                           sapply(multihitPositions, length)),
          "sampleID"=metadata[j,"sampleName"],
          "position"=start(flank(unlist(multihitPositions), width=-1, start=TRUE, both=FALSE)),
          "chr"=as.character(seqnames(unlist(multihitPositions))),
          "strand"=as.character(strand(unlist(multihitPositions))) )
        
        # change to the right class as database
        multihitPositions$multihitID <- as(multihitPositions$multihitID, "integer")
        multihitPositions$sampleID <- as(multihitPositions$sampleID, "integer")
        multihitPositions$position <- as(multihitPositions$position, "integer")
        multihitPositions$chr <- as(multihitPositions$chr, "character")
        multihitPositions$strand <- as(multihitPositions$strand, "character")
        
        multihitLengths <- data.frame(
          "multihitID"=rep(seq(length(multihitLengths))+currentMaxMultihitID,
                           sapply(multihitLengths, nrow)),
          "length"=as.integer(as.character(do.call(rbind, multihitLengths)$Var1)),
          "count"=do.call(rbind, multihitLengths)$Freq )
        
        # change to the right class as database
        multihitLengths$multihitID <- as(multihitLengths$multihitID, "integer")
        multihitLengths$length <- as(multihitLengths$length, "integer")
        multihitLengths$count <- as(multihitLengths$count, "integer")
        
        # MINE: instead of loading it into the database, Ill read them here. 
        # assign(paste0("multihitpositions_for_",file), multihitPositions ,envir = .GlobalEnv)
        # assign(paste0("multihitlengths_for_",file), multihitLengths ,envir = .GlobalEnv)
        # 
        # load table multihitpositions
        message("Loading multihitpositions:", nrow(multihitPositions), " entries")
        ### stopifnot( dbWriteTable(dbConn, "multihitpositions", multihitPositions, append=T, row.names=F) )
        #null <- apply(multihitPositions, 1, insertSQL, dbConn=dbConn, colNames=colnames(multihitPositions), table='multihitPositions')            
        
        # load table multihitlengths
        message("Loading multihitlengths: ", nrow(multihitLengths), " entries")
        ### stopifnot( dbWriteTable(dbConn, "multihitlengths", multihitLengths, append=T, row.names=F) )
        #null <- apply(multihitLengths, 1, insertSQL, dbConn=dbConn, colNames=colnames(multihitLengths), table='multihitLengths')
        
        newMaxMultihitID = currentMaxMultihitID + length(unique(multihitPositions$multihitID))
        currentMaxMultihitID <- newMaxMultihitID  
      }
    }
    
  }
  

}



################################################################################
save.image (file =file.path(RESULTS_DIR,"Inspiired_CD33_jun.RData"))
savehistory(file =file.path(RESULTS_DIR,"Inspiired_CD33_jun.Rhistory"))
sink(file =file.path(RESULTS_DIR,"Inspiired_CD33_jun.txt"))
toLatex(sessionInfo())
sink(NULL)
################################################################################
