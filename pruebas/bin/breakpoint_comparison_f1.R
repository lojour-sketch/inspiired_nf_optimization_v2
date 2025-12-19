#!/usr/bin/env Rscript

library(GenomicRanges)
library(regioneR)
library(GenomeInfoDb)
library(readr)
library(VennDiagram)
library(openxlsx)

# INPUT FILES
samples <- list(
  "HD68_KO" = list(
    insp = "/home/lrenteria/inspiired_nf/maren_inspiired_run/CD33_worked/CD33_HD68_march/HD68_KO/allSites.RData",
    nf   = "/home/lrenteria/inspiired_nf/results/17_sitesfinal_inspiired/CD33/HD68_KO_CART_S2_allsites_nostandard.rds"
  ),
  "HD68_wt" = list(
    insp = "/home/lrenteria/inspiired_nf/maren_inspiired_run/CD33_worked/CD33_HD68_march/HD68_wt/allSites.RData",
    nf   = "/home/lrenteria/inspiired_nf/results/17_sitesfinal_inspiired/CD33/HD68_wt_CART_S1_allsites_nostandard.rds"
  ),
  "HD70_KO" = list(
    insp = "/home/lrenteria/inspiired_nf/maren_inspiired_run/CD33_worked/CD33_HD70_march/HD70_KO/allSites.RData",
    nf   = "/home/lrenteria/inspiired_nf/results/17_sitesfinal_inspiired/CD33/HD70_KO_CART_S4_allsites_nostandard.rds"
  ),
  "HD70_wt" = list(
    insp = "/home/lrenteria/inspiired_nf/maren_inspiired_run/CD33_worked/CD33_HD70_march/HD70_wt/allSites.RData",
    nf   = "/home/lrenteria/inspiired_nf/results/17_sitesfinal_inspiired/CD33/HD70_wt_CART_S3_allsites_nostandard.rds"
  )
)

# PARAMETERS
tolerance <- 10
perm_num <- 100

# HELPER FUNCTIONS
load_granges <- function(path, type="rds") {
  if(type == "rds") {
    readRDS(path)
  } else {
    preload <- load(path)
    get(preload)
  }
}

breakpoints <- function(gr) {
  pos <- ifelse(strand(gr) == "+", start(gr), end(gr))
  GRanges(seqnames = seqnames(gr), ranges = IRanges(pos, pos))
}

compute_metrics <- function(insp_gr, nf_gr, tolerance) {
  insp_bp <- breakpoints(insp_gr)
  nf_bp   <- breakpoints(nf_gr)

  insp_exp <- resize(insp_bp, width = 2*tolerance+1, fix="center")
  nf_exp <- resize(nf_bp, width = 2*tolerance+1, fix="center")

  recall <- sum(countOverlaps(insp_exp, nf_bp) > 0) / length(insp_bp)
  precision <- sum(countOverlaps(nf_exp, insp_bp) > 0) / length(nf_bp)
  F1 <- 2 * (precision * recall) / (precision + recall)

  list(
    recall = recall,
    precision = precision,
    F1 = F1,
    insp_total = length(insp_bp),
    nf_total = length(nf_bp),
    insp_dup = sum(duplicated(insp_bp)),
    nf_dup = sum(duplicated(nf_bp)),
    insp_unique = length(unique(insp_bp)),
    nf_unique = length(unique(nf_bp))
  )
}

# RUN ANALYSIS FOR ALL SAMPLES
all_metrics <- list()
for(samp_name in names(samples)) {
  cat("Processing:", samp_name, "\n")
  insp <- load_granges(samples[[samp_name]]$insp, type="RData")
  nf   <- load_granges(samples[[samp_name]]$nf, type="rds")

  metrics <- compute_metrics(insp, nf, tolerance)
  metrics$sample <- samp_name
  all_metrics[[samp_name]] <- metrics

  # Venn Diagram
  insp_bp <- unique(breakpoints(insp))
  nf_bp <- unique(breakpoints(nf))
  venn_file <- paste0(samp_name, "_venn.pdf")
  pdf(venn_file)
  draw.pairwise.venn(area1 = length(insp_bp), area2 = length(nf_bp), 
                     cross.area = sum(countOverlaps(insp_bp, nf_bp) > 0),
                     category = c("Inspiired", "NF"), 
                     fill = c("blue", "red"))
  dev.off()
}

# EXPORT RESULTS TO EXCEL
metrics_df <- do.call(rbind, lapply(all_metrics, as.data.frame))
write.xlsx(metrics_df, file="comparison_metrics.xlsx")
message("Analysis complete. Metrics saved in comparison_metrics.xlsx")
