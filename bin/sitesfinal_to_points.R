#!/usr/bin/env Rscript

# ------------------------------------------------------------------
# Author: Libe Renteria Aizpurua
# Date: 2026-01-07 
#
# This script annotates insertion sites represented as points using ChIPseeker and generates plots and an Excel file with the annotations.
#
# ------------------------------------------------------------------


# LOAD REQUIRED LIBRARIES
library(GenomicRanges)
library(ChIPseeker)
library(WriteXLS)
library(rtracklayer)
library(GenomicFeatures)
library(GenomeInfoDb)
library(dplyr)
library(org.Hs.eg.db)
library(clusterProfiler)

#debug
message("=== STARTING ANNOTATION SCRIPT WITH POINT INFORMATION ===")
message("Timestamp: ", Sys.time())

# Get command line arguments
args <- commandArgs(trailingOnly = TRUE)
sample <- args[1] 
sitesFinal <- args[2]
ref_genome <- args[3]
txdbFile <- args[4]

# Force writable cache locations inside task work dir (containers may have read-only /root)
cache_root <- file.path(getwd(), ".clusterprofiler_cache")
dir.create(cache_root, recursive = TRUE, showWarnings = FALSE)
xdg_data_home <- file.path(cache_root, "xdg_data")
xdg_cache_home <- file.path(cache_root, "xdg_cache")
dir.create(xdg_data_home, recursive = TRUE, showWarnings = FALSE)
dir.create(xdg_cache_home, recursive = TRUE, showWarnings = FALSE)
Sys.setenv(HOME = cache_root, XDG_DATA_HOME = xdg_data_home, XDG_CACHE_HOME = xdg_cache_home)
message("Using writable local cache at: ", xdg_data_home)

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
} else if (txdbFile == "TxDb.Hsapiens.UCSC.hg18.knownGene") {
    library(TxDb.Hsapiens.UCSC.hg18.knownGene)
    txdb <- TxDb.Hsapiens.UCSC.hg18.knownGene
    message("Loaded hg18 knownGene")
} else {
    stop("Unsupported TxDb file: ", txdbFile)
}

if(endsWith(sitesFinal, ".RData")){
    message("  Reading sitesFinal: ", sitesFinal)
    object <- load(sitesFinal)
    sites.final <- get(object)
} else {
    message("  Reading sitesFinal: ", sitesFinal)
    sites.final <- readRDS(sitesFinal)
    message("  sites.final dimensions: ", length(sites.final))
}

#stop the script if sitesfinal is just one line
if(length(sites.final) == 1){
    message("Sample has no sites to annotate")
    
    # Crea outputs placeholder
    pdf(paste0("fig_points_", sample, ".pdf"), width = 8, height = 6)
    plot.new()
    text(0.5, 0.5, 
         paste0("Sample: ", sample, "\n\nNo sites to annotate"),
         cex = 1.2)
    dev.off()
    
    WriteXLS::WriteXLS(
        x = list(
            Summary = data.frame(
                Sample = sample,
                Status = "No annotation",
                Reason = "sites.final had no sites"
            )
        ),
        ExcelFileName = paste0("annotated_points_", sample, ".xlsx"),
        SheetNames = "Summary"
    )
    
    # Sale con éxito (exit code 0)
    quit(save = "no", status = 0)
}

#sites.final is a GRanges object, so we convert it to a data.frame first
df <- as.data.frame(mcols(sites.final))
df$seqnames <- as.character(seqnames(sites.final))
df$start <- start(sites.final)
df$end <- end(sites.final)
df$strand <- as.character(strand(sites.final))


#now we will collapse the rows that have the same insertion in the same strand
#for insertions on the + strand, we will collapse those that have the same start (insertion of r2)
#for insertions on the - strand, we will collapse those that have the same end (insertion of r2)
df_collapsed <- df %>%
  mutate(r2_pos = if_else(strand == "+", start, end)) %>%
  group_by(seqnames, strand, r2_pos) %>%
  summarise(
    counts = as.integer(n()),
    strand = first(strand),
    start = first(start),
    end = first(end),
    revmap = paste(sort(unique(unlist(revmap))), collapse=";"),
    pairingID = first(pairingID),
    samplename = first(samplename),
  ) %>%
  ungroup()

# we will reconstruct a GRanges object to be able to plot with ChipSeek
sitescollapsed <- GRanges(
  seqnames = df_collapsed$seqnames,
  ranges   = IRanges(start=df_collapsed$start, end=df_collapsed$end),
  strand   = df_collapsed$strand,
  counts     = df_collapsed$counts,
  revmap     = df_collapsed$revmap,
  pairingID  = df_collapsed$pairingID,
  samplename = df_collapsed$samplename
)

#annotation and plotting
peak <- sitescollapsed
    chromosomes = c(paste0("chr",1:22), "chrX", "chrY")

    message("Debugging counts column:")

        stopifnot("counts" %in% names(mcols(peak)))

        mcols(peak)$counts <- as.numeric(mcols(peak)$counts)
        message("counts column class: ", class(mcols(peak)$counts))
        message("counts column length: ", length(mcols(peak)$counts))

        print(mcols(peak))
        print(names(mcols(peak)))

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

    pdf(paste0("fig_points_", sample, '.pdf'))
    print(covplot(peak, weightCol = "counts", chrs = chromosomes,
                title = paste0("Insertion Sites over Chromosomes sample: ", sample)))
    print(plotAnnoPie(peakAnno))
    print(plotAnnoBar(peakAnno))
    print(upsetplot(peakAnno))
    print(vennpie(peakAnno))
    print(upsetplot(peakAnno, vennpie = TRUE))
    # Only plot if there are both upstream and downstream distances
    distances <- as.data.frame(peakAnno)$distanceToTSS
    if (sum(!is.na(distances)) > 0 && length(unique(sign(distances))) > 1) {
        print(plotDistToTSS(peakAnno))
    } else {
        message("Skipping plotDistToTSS: not enough data on both sides of TSS")
    }
    peakAnno.dfr <- as.data.frame(peakAnno)
    entrez_ids <- peakAnno.dfr$geneId
    entrez_ids <- entrez_ids[!is.na(entrez_ids)]
    entrez_ids <- as.character(entrez_ids)

    if(length(entrez_ids) > 0) {
        message("Calculating GO and KEGG enrichment...")
        goenrichment <- enrichGO(
            gene = entrez_ids,
            OrgDb = org.Hs.eg.db,
            keyType = "ENTREZID",
            ont = "BP",
            pAdjustMethod = "BH",
            qvalueCutoff = 0.05,
            readable = TRUE
        )

        if(!is.null(goenrichment) && nrow(as.data.frame(goenrichment)) > 0) {
            message("GO enrichment results:")
            print(dotplot(goenrichment, showCategory = 50, font.size = 4))
        } else {
            message("No GO enrichment results found, skipping dotplot.")
        }

        keggenrichment <- tryCatch(
            enrichKEGG(
                gene = entrez_ids,
                organism = "hsa",
                pvalueCutoff = 0.05
            ),
            error = function(e) {
                message("KEGG enrichment failed; continuing without KEGG plot. Error: ", conditionMessage(e))
                NULL
            }
        )

        if(!is.null(keggenrichment) && nrow(as.data.frame(keggenrichment)) > 0) {
            message("KEGG enrichment results:")
            print(dotplot(keggenrichment))
        } else {
            message("No KEGG enrichment results found, skipping dotplot.")
        }
    } else {
        message("No ENTREZ IDs found for enrichment analysis, skipping GO and KEGG enrichment.")
    }

    message("PDF generation completed")
    dev.off()


# Add the missing columns from sites.final
final_output <- peakAnno.dfr

#save the complete Excel
WriteXLS::WriteXLS(final_output, paste0("annotated_points_", sample, '.xlsx'))
