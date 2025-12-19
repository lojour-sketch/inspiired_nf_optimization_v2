#!/usr/bin/env Rscript

#load libraries
library(readxl)
library(GenomicRanges)
library(VennDiagram)

# Get command line arguments
args <- commandArgs(trailingOnly = TRUE)
annot_nf <- args[1]
annot_inspiired <- args[2]


#we read the annotation xlsx files
nf_annot <- read_excel(annot_nf)
insp_annot <- read_excel(annot_inspiired)

#convert the files to granges
nf_annot_gr <- GRanges(seqnames = nf_annot$seqnames,
                          ranges = IRanges(start = nf_annot$start, end = nf_annot$end),
                          strand = nf_annot$strand)
insp_annot_gr <- GRanges(seqnames = insp_annot$seqnames,
                          ranges = IRanges(start = insp_annot$start, end = insp_annot$end),
                          strand = insp_annot$strand)

#now we find the overlaping sites
overlaping_sites <- findOverlaps(nf_annot_gr, insp_annot_gr)
shared <- nf_annot[queryHits(overlaping_sites), ]
unique_nf <- nf_annot[-unique(queryHits(overlaping_sites)), ]
unique_insp <- insp_annot[-unique(subjectHits(overlaping_sites)), ]

#create counts for venn diagram
venn_counts <- c(
    nf = nrow(nf_annot),
    insp = nrow(insp_annot),
    shared = length(unique(queryHits(overlaping_sites)))
)

# Plot Venn diagram
venn.plot <- draw.pairwise.venn(
  area1 = venn_counts["nf"],
  area2 = venn_counts["insp"],
  cross.area = venn_counts["shared"],
  category = c("NF", "INSPIIRED"),
  fill = c("skyblue", "orange"),
  alpha = 0.5,
  cat.cex = 1.2
)
