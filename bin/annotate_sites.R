#!/usr/bin/env Rscript

# LOAD REQUIRED LIBRARIES
library(GenomicRanges)
library(ChIPseeker)
library(WriteXLS)
library(rtracklayer)
library(GenomicFeatures)
library(GenomeInfoDb)
library(dplyr)
library(org.Hs.eg.db)

#debug
message("=== STARTING ANNOTATION SCRIPT ===")
message("Timestamp: ", Sys.time())

# Get command line arguments
args <- commandArgs(trailingOnly = TRUE)
sample <- args[1]
sitesFinal <- args[2] 
allSites <- args[3]
ref_genome <- args[4]
txdbFile <- args[5]

#debug
message("Command line arguments:")
message("  sample: ", sample)
message("  sitesFinal: ", sitesFinal)
message("  allSites: ", allSites)
message("  ref_genome: ", ref_genome)
message("  txdbFile: ", txdbFile)

# Load the TxDb package correctly
message("Loading TxDb package: ", txdbFile)
if (txdbFile == "TxDb.Hsapiens.UCSC.hg19.knownGene") {
    library(TxDb.Hsapiens.UCSC.hg19.knownGene)
    txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
    message("Loaded hg19 knownGene")
} else if (txdbFile == "TxDb.Hsapiens.UCSC.hg38.refGene") {
    library(TxDb.Hsapiens.UCSC.hg38.refGene)
    txdb <- TxDb.Hsapiens.UCSC.hg38.refGene
    message("Loaded hg38 refGene")
} else {
    stop("Unsupported TxDb file: ", txdbFile)
}

message("Reading RDS files...")
message("  Reading allSites: ", allSites)
allsites_gr <- readRDS(allSites)
message("  allsites_gr dimensions: ", length(allsites_gr))

message("  Reading sitesFinal: ", sitesFinal)
sites.final <- readRDS(sitesFinal)
message("  sites.final dimensions: ", length(sites.final))


#we will create a column with the complete chromosome name and another with just the number
mcols(sites.final)$orig_chr <- as.character(seqnames(sites.final))
seqlevelsStyle(sites.final) <- "UCSC"

if(length(sites.final) > 0) {
    message("Processing ", length(sites.final), " sites...")

    #we will create sites dataframe with condensed sites
    message("Creating sites dataframe...")
    sites <- data.frame(
        "siteID" = seq(length(sites.final)),
        "sampleID" = sample,
        "position" = start(flank(sites.final, -1, start = T)),
        "chr" = as.character(seqnames(sites.final)),
        "strand" = as.character(strand(sites.final))
    )

    ## change to the right class as database
    sites$siteID <- as(sites$siteID, "integer")
    sites$sampleID <- as(sites$sampleID, "character")
    sites$position <- as(sites$position, "integer")
    sites$chr <- as(sites$chr, "character")
    sites$strand <- as(sites$strand, "character")

    # Reorder allSites to match sites.final
    message("Reordering allSites...")
    allSites <- allsites_gr[unlist(sites.final$revmap)]
    message("  allSites after reordering: ", length(allSites))
    #as we saw that some rows are repeated, we will make the granges unique
    #we think this is caused because in the standardization step some sites are grouped together that then are not merged in the dereplication step
    message("Making sites unique...")
    allSites <- unique(allSites)
    sites.final <- unique(sites.final)
    message("  Unique allSites: ", length(allSites))
    message("  Unique sites.final: ", length(sites.final))

    message("Calculating PCR breakpoints...")
    pcrBreakpoints <- sort(paste0(as.integer(Rle(sites$siteID, sapply(sites.final$revmap, length))),
                            ".",
                            start(flank(allSites, -1, start = F))))

    condensedPCRBreakpoints <- strsplit(unique(pcrBreakpoints), "\\.")
    #with pcrbreakpoints we save all the insertion sites that were clustered together
    pcrBreakpoints <- data.frame(
        "siteID" = sapply(condensedPCRBreakpoints, "[[", 1),
        "breakpoint" = sapply(condensedPCRBreakpoints, "[[", 2),
        "count" = runLength(Rle(match(pcrBreakpoints, unique(pcrBreakpoints))))
    )
    message("Created Sites and PCRbreakpoints objects")
    

    # MAREN: instead of loading it into the database, Ill read them here. 
    assign(paste0("Sites_for_",sample), sites ,envir = .GlobalEnv)
    assign(paste0("PCRbreakpoints_for_",sample), pcrBreakpoints ,envir = .GlobalEnv)
}

#plotting
if(length(sites.final) > 0) {
    peak <- sites.final
    chromosomes = c(paste0("chr",1:22), "chrX", "chrY")

    message("Starting peak annotation...")
    message("  Number of peaks to annotate: ", length(peak))
    message("  Peak object size: ", format(object.size(peak), units = "MB"))


    start_time <- Sys.time()
    peakAnno <- annotatePeak(peak, tssRegion = c(-3000, 3000),
                            TxDb = txdb, annoDb="org.Hs.eg.db")
    end_time <- Sys.time()

    annotation_time <- end_time - start_time
    message("Peak annotation completed in ", round(annotation_time, 2), " ", units(annotation_time))

    message("Generating plots...")

    pdf(paste0("fig_annot_", sample, '.pdf'))
    print(covplot(peak, weightCol = "counts", chrs = chromosomes,
                title = paste0("Insertion Sites over Chromosomes sample: ", sample)))
    print(plotAnnoPie(peakAnno))
    print(plotAnnoBar(peakAnno))
    print(upsetplot(peakAnno))
    print(vennpie(peakAnno))
    print(upsetplot(peakAnno, vennpie = TRUE))
    print(plotDistToTSS(peakAnno))
    message("PDF generation completed")
    dev.off()
}
# to create exactly the file with the columns we want
message("Creating output data frame...")
peakAnno.dfr <- as.data.frame(peakAnno)

#debugging:
message("Peak annotation data frame structure:")
str(peakAnno.dfr)
message("Peak annotation data frame head:")
head(peakAnno.dfr)
message("Peak annotation data frame columns:")
colnames(peakAnno.dfr)


# Añadir las columnas que faltan desde sites.final
final_output <- peakAnno.dfr


  
  # Verificar que tenemos todas las columnas
message("Columnas finales: ", paste(names(final_output), collapse = ", "))
message("Número de columnas: ", ncol(final_output))

# Guardar el Excel completo
WriteXLS::WriteXLS(final_output, paste0("annotated_", sample, '.xlsx'))



message("Analysis complete for sample: ", sample)