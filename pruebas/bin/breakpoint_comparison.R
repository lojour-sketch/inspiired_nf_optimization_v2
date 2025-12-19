#!/usr/bin/env Rscript

# LOAD REQUIRED LIBRARIES
library(GenomicRanges)
library(regioneR)
library(GenomeInfoDb)
library(readr)

################################# INPUTS #################################
        HD68KO_allSites_nf <- "/home/lrenteria/inspiired_nf/results/17_sitesfinal_inspiired/CD33/HD68_KO_CART_S2_allsites_nostandard.rds"
        HD68KO_allSites_insp <- "/home/lrenteria/inspiired_nf/maren_inspiired_run/CD33_worked/CD33_HD68_march/HD68_KO/allSites.RData"

        HD68wt_allSites_nf <- "/home/lrenteria/inspiired_nf/results/17_sitesfinal_inspiired/CD33/HD68_wt_CART_S1_allsites_nostandard.rds"
        HD68wt_allSites_insp <- "/home/lrenteria/inspiired_nf/maren_inspiired_run/CD33_worked/CD33_HD68_march/HD68_wt/allSites.RData"

        HD70KO_allSites_nf <- "/home/lrenteria/inspiired_nf/results/17_sitesfinal_inspiired/CD33/HD70_KO_CART_S4_allsites_nostandard.rds"
        HD70KO_allSites_insp <- "/home/lrenteria/inspiired_nf/maren_inspiired_run/CD33_worked/CD33_HD70_march/HD70_KO/allSites.RData"

        HD70wt_allSites_nf <- "/home/lrenteria/inspiired_nf/results/17_sitesfinal_inspiired/CD33/HD70_wt_CART_S3_allsites_nostandard.rds"
        HD70wt_allSites_insp <- "/home/lrenteria/inspiired_nf/maren_inspiired_run/CD33_worked/CD33_HD70_march/HD70_wt/allSites.RData"


        #load allSites
        HD68_KO_allSites_nf <- readRDS(HD68KO_allSites_nf)
        preload <- load(HD68KO_allSites_insp)
        HD68_KO_allSites_insp <- get(preload)

        #load sites.final
        HD68_wt_allSites_nf <- readRDS(HD68wt_allSites_nf)
        preload <- load(HD68wt_allSites_insp)
        HD68_wt_allSites_insp <- get(preload)

        #load sites.final
        HD70_KO_allSites_nf <- readRDS(HD70KO_allSites_nf)
        preload <- load(HD70KO_allSites_insp)
        HD70_KO_allSites_insp <- get(preload)


        #load sites.final
        HD70_wt_allSites_nf <- readRDS(HD70wt_allSites_nf)
        preload <- load(HD70wt_allSites_insp)
        HD70_wt_allSites_insp <- get(preload)

################################# PARAMETERS #################################
tolerance <- 5
perm_num <- 1000

breakpoints <- function(gr) {
  pos <- ifelse(strand(gr) == "+", start(gr), end(gr))
  GRanges(seqnames = seqnames(gr), ranges = IRanges(pos, pos))
}

compare_datasets <- function(insp, nf, sample_name) {
  insp_bp <- breakpoints(insp)
  nf_bp   <- breakpoints(nf)

  insp_exp <- resize(insp_bp, width = 2 * tolerance + 1, fix = "center")

  hits <- sum(countOverlaps(insp_exp, nf_bp) > 0)
  insp_proportion <- hits / length(insp_bp)

  ref <- insp_bp
  eval_fun <- function(x) {
    sum(countOverlaps(resize(x, 2 * tolerance + 1, "center"), ref) > 0) / length(x)
  }

  pt <- permTest(
    A = nf_bp,
    ntimes = perm_num,
    randomize.function = randomizeRegions,
    evaluate.function = eval_fun,
    alternative = "greater"
  )

  message("\n========== ", sample_name, " ==========")
  message("Proportion of reference breakpoints recovered: ", insp_proportion)
  message("Empirical p-value from permutation test: ", pt$eval_fun$pval)
  # Recall: fraction of reference sites recovered
        recall <- hits / length(insp_bp)
        message("Recall: ", recall)
  # Precision: fraction of dataset sites that match reference
        nf_exp <- resize(nf_bp, width = 2 * tolerance + 1, fix = "center")
        hits_precision <- sum(countOverlaps(nf_exp, insp_bp) > 0)
        precision <- hits_precision / length(nf_bp)
# F1-score: harmonic mean
        F1 <- 2 * (precision * recall) / (precision + recall)
        message("F1-score: ", F1)

}

################################# RUN FOR ALL SAMPLES #################################
compare_datasets(HD68_KO_allSites_insp, HD68_KO_allSites_nf, "HD68_KO")
compare_datasets(HD68_wt_allSites_insp, HD68_wt_allSites_nf, "HD68_wt")
compare_datasets(HD70_KO_allSites_insp, HD70_KO_allSites_nf, "HD70_KO")
compare_datasets(HD70_wt_allSites_insp, HD70_wt_allSites_nf, "HD70_wt")