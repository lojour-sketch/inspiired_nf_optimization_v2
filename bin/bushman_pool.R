#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly=TRUE)
if(!length(args) %in% c(7L,8L)) stop("Usage: SAMPLES CALLS_DIRECTORY SOURCES WINDOW METHOD ANNOTATION GENOME [FAILURE_POLICY]")
failure_policy <- if(length(args)==8L) args[8] else "fail"
if(!failure_policy %in% c("fail","record")) stop("Invalid abundance failure policy")
own <- dirname(normalizePath(sub("^--file=","",grep("^--file=",commandArgs(),value=TRUE)[1])))
source(file.path(own,"bushman_runtime.R"))
e <- bushman_runtime(args[3],as.integer(args[4]))
metadata <- data.table::fread(args[1],colClasses="character")
scope <- strsplit(args[7],"__",fixed=TRUE)[[1]]
genome_name <- scope[1]; report_group <- paste(scope[-1],collapse="__")
metadata <- metadata[refGenome==genome_name & Report_Group==report_group]
method <- args[5]; annotation <- args[6]
if(!method %in% c("sonic","fragment_diversity")) stop("Invalid abundance method")
source(file.path(args[3],"geneTherapyPatientReportMaker","estimatedAbundance.R"),local=e)
source(file.path(args[3],"geneTherapyPatientReportMaker","populationInfo.R"),local=e)
started <- proc.time()
all <- list(); multi <- list(); summaries <- list()
for(i in seq_len(nrow(metadata))) {
    row <- metadata[i]; path <- file.path(args[2],row$Sample_ID)
    if(!file.exists(file.path(path,"allSites.RData"))) stop("Missing expected sample: ",row$Sample_ID)
    sites <- get(load(file.path(path,"allSites.RData")))
    if(length(sites)) {
        # Match the DB's unique PCR-breakpoint rows per library, not read counts.
        sites <- sites[!duplicated(as.character(granges(sites)))]
        mcols(sites) <- S4Vectors::DataFrame(sampleName=row$Sample_ID,GTSP=row$Specimen_ID,
            Replicate_ID=row$Replicate_ID)
        all[[row$Sample_ID]] <- sites
    }
    mh <- get(load(file.path(path,"multihitData.RData")))
    for(j in seq_along(mh$clusteredMultihitLengths)) {
        lens <- mh$clusteredMultihitLengths[[j]]
        multi[[paste(row$Sample_ID,j,sep=":")]] <- data.frame(sampleName=row$Sample_ID,
            GTSP=row$Specimen_ID,multihitID=paste(row$Sample_ID,j,sep=":"),
            length=as.integer(as.character(lens$Var1)),estAbund=nrow(lens),
            method="original_multihit_distinct_lengths")
    }
    summaries[[row$Sample_ID]] <- as.data.frame(get(load(file.path(path,"stats.RData"))))
}
unique_sites <- if(length(all)) do.call(c,unname(all)) else GenomicRanges::GRanges()
raw <- unique_sites
if(length(unique_sites)) unique_sites <- e$standardizeSites(unique_sites)
if(length(unique_sites)) unique_sites$posid <- paste0(seqnames(unique_sites),strand(unique_sites),
    start(flank(unique_sites,-1,start=TRUE)))
if(length(unique_sites)) {
    unique_sites$replicate <- 0L
    for(specimen in unique(unique_sites$GTSP)) {
        selected <- which(unique_sites$GTSP==specimen)
        unique_sites$replicate[selected] <- as.integer(as.factor(unique_sites$sampleName[selected]))
    }
}
saveRDS(list(raw=raw,standardized=unique_sites,window_bp=as.integer(args[4]),
    scope="all selected specimens in report group/reference before specimen splitting"),"pooled_sites.rds")
estimates <- list(); fits <- list(); statuses <- list(); occupancy <- list()
e$estAbund <- function(...) {
    fit <- sonicLength::estAbund(...)
    fits[[specimen]] <<- fit
    fit
}
for(specimen in unique(metadata$Specimen_ID)) {
    sites <- unique_sites[unique_sites$GTSP==specimen]
    if(!length(sites)) {statuses[[specimen]] <- data.frame(GTSP=specimen,status="empty",warning="");next}
    # Original report derives technical strata from each sampleName.
    sites$replicate <- as.integer(as.factor(sites$sampleName))
    locations <- paste0(seqnames(sites),strand(sites),start(flank(sites,-1,start=TRUE)))
    input <- data.frame(GTSP=specimen,location=locations,fragment_length=width(sites),replicate=sites$replicate)
    occupancy[[specimen]] <- input
    warning_text <- character()
    result <- tryCatch(withCallingHandlers(
        e$getEstimatedAbundance(sites,use.sonicLength=(method=="fragment_diversity")),
        warning=function(w){warning_text <<- c(warning_text,conditionMessage(w));invokeRestart("muffleWarning")}),
        error=function(err){warning_text <<- c(warning_text,conditionMessage(err));NULL})
    if(is.null(result)) {
        statuses[[specimen]] <- data.frame(GTSP=specimen,status="estimator_failed",warning=paste(warning_text,collapse="; "))
        next
    }
    # Archive a direct fit only for model inspection; original rounded reporting
    # above determines abundance, proportion and rank. Multihits never enter it.
    result$GTSP <- specimen; result$abundance_method <- method
    estimates[[specimen]] <- result
    statuses[[specimen]] <- data.frame(GTSP=specimen,status=if(length(warning_text)) "estimated_with_warnings" else "estimated",
        warning=paste(warning_text,collapse="; "))
}
final <- if(length(estimates)) do.call(c,unname(estimates)) else GenomicRanges::GRanges()
failed_specimens <- names(statuses)[vapply(statuses,function(x) x$status=="estimator_failed",logical(1))]
report_state <- if(length(failed_specimens)) "completed_with_unavailable_abundance" else "completed"
if(length(final) && "Timepoint" %in% names(metadata) &&
    !(failure_policy=="record" && length(failed_specimens))) {
    clinical <- c("Timepoint","CellType","Patient","Trial","VCN")
    for(field in clinical[clinical %in% names(metadata)]) {
        by_specimen <- split(metadata[[field]],metadata$Specimen_ID)
        if(any(vapply(by_specimen,function(x) length(unique(x))!=1L,logical(1))))
            stop("Conflicting specimen metadata: ",field)
        values <- vapply(by_specimen,function(x) x[1],character(1))
        mcols(final)[[field]] <- values[final$GTSP]
        mcols(unique_sites)[[field]] <- values[unique_sites$GTSP]
    }
    e$diversity <- vegan::diversity; e$estimateR <- vegan::estimateR
    e$acast <- reshape2::acast; e$gini <- reldist::gini
    set.seed(1L) # Record the seed for the original jackknife's random assignment.
    population <- e$getPopulationInfo(unique_sites,final,"GTSP")
    population$Replicates <- vapply(split(unique_sites$replicate,unique_sites$GTSP),max,integer(1))
    population$UniqueSites <- vapply(split(final$posid,final$GTSP),function(x) length(unique(x)),integer(1))
    population$InferredCells <- vapply(split(unique_sites,unique_sites$GTSP),function(x)
        nrow(unique(data.frame(replicate=x$replicate,posid=x$posid,width=width(x)))),integer(1))
    data.table::fwrite(population,"original_population_by_specimen.tsv",sep="\t")
    timepoints <- e$getPopulationInfo(unique_sites,final,"Timepoint")
    timepoints$UniqueSites <- vapply(split(final$posid,final$Timepoint),function(x) length(unique(x)),integer(1))
    data.table::fwrite(timepoints,"original_population_by_timepoint.tsv",sep="\t")
    data.table::fwrite(data.frame(status="computed",method="frozen Bushman populationInfo.R",jackknife_seed=1L),
        "original_population_status.tsv",sep="\t")
} else data.table::fwrite(data.frame(status="unavailable",reason=if(failure_policy=="record" && length(failed_specimens))
        "Population summaries unavailable: at least one specimen has no valid abundance estimate" else
        "Timepoint metadata or estimated unique sites absent"),
    "original_population_status.tsv",sep="\t")
saveRDS(final,"original_unique_abundance.rds")
saveRDS(fits,"original_sonic_fits.rds")
if(length(final)) {
    points <- flank(final,-1,start=TRUE)
    if(annotation=="bushman") {
        for(genome in unique(metadata$refGenome)) {
            ref <- unique(metadata[refGenome==genome]$refseq_rds)
            if(length(ref)!=1L) stop("One frozen RefSeq annotation per assembly")
            genes <- readRDS(ref)
            points <- hiAnnotator::getNearestFeature(points,genes,colnam="nearest_refSeq_gene",feature.colnam="name2")
            points <- hiAnnotator::getNearestFeature(points,genes,colnam="nearest_refSeq_gene",side="5p",feature.colnam="name2")
            points <- hiAnnotator::getSitesInFeature(points,genes,colnam="inGene",feature.colnam="name2")
            oncofile <- if(grepl("^mm",genome)) "allonco_no_pipes.mm.csv" else "allonco_no_pipes.csv"
            onco <- scan(file.path(args[3],"geneTherapyPatientReportMaker",oncofile),what="character",quiet=TRUE)
            onco <- onco[!grepl("geneName",onco,ignore.case=TRUE)]
            points <- hiAnnotator::getNearestFeature(points,genes[toupper(genes$name2)%in%toupper(onco)],
                colnam="NrstOnco",side="5p",feature.colnam="name2")
            adverse <- read.delim(file.path(args[3],"geneTherapyPatientReportMaker","humanLymph.tsv"))$symbol
            points$geneMark <- ""
            in_adverse <- vapply(strsplit(points$inGene,","),function(x) any(x%in%adverse),logical(1))
            in_onco <- vapply(strsplit(points$inGene,","),function(x) any(x%in%onco),logical(1))
            points$geneMark <- ifelse(points$nearest_refSeq_gene%in%adverse|in_adverse,"!","")
            near_onco <- points$nearest_refSeq_gene%in%onco & abs(points$nearest_refSeq_geneDist)<50000
            points$geneMark <- paste0(points$geneMark,ifelse(near_onco|in_onco,"~",""),
                ifelse(toupper(points$inGene)!="FALSE","*",""))
            points$nearest_refSeq_geneUnmarked <- points$nearest_refSeq_gene
            points$nearest_refSeq_gene <- paste0(points$nearest_refSeq_gene,points$geneMark)
        }
    }
    saveRDS(points,"original_annotated_sites.rds")
    # Specimen-local GRanges labels can repeat after pooling. They are not TSV
    # columns; use automatic table row names and preserve the scientific object.
    data.table::fwrite(as.data.frame(points,row.names=NULL),"original_unique_abundance.tsv",sep="\t")
} else data.table::fwrite(data.frame(GTSP=character(),posid=character(),estAbund=numeric(),estAbundProp=numeric()),
    "original_unique_abundance.tsv",sep="\t")
multi_table <- if(length(multi)) data.table::rbindlist(multi) else
    data.table(sampleName=character(),GTSP=character(),multihitID=character(),length=integer(),estAbund=integer(),method=character())
if(nrow(multi_table)) {
    multi_table[,replicate:=as.integer(as.factor(sampleName)),by=multihitID]
    for(field in intersect(c("Patient","Timepoint","CellType"),names(metadata)))
        multi_table[[field]] <- metadata[[field]][match(multi_table$sampleName,metadata$Sample_ID)]
    if(all(c("Patient","Timepoint","CellType") %in% names(multi_table)))
        multi_table[,Rank:=rank(-estAbund,ties.method="max"),by=.(Patient,Timepoint,CellType)]
}
data.table::fwrite(multi_table,"original_multihit_abundance.tsv",sep="\t")
data.table::fwrite(data.table::rbindlist(statuses),"abundance_status.tsv",sep="\t")
data.table::fwrite(data.table::rbindlist(occupancy),"original_estimator_input.tsv",sep="\t")
data.table::fwrite(data.table::rbindlist(summaries,fill=TRUE),"sample_stats.tsv",sep="\t")
jsonlite::write_json(list(elapsed_seconds=unname((proc.time()-started)[3]),window_bp=as.integer(args[4]),
    abundance_method=method,abundance_failure_policy=failure_policy,report_state=report_state,
    unavailable_specimens=unname(failed_specimens),normalization="rounded unique-site theta within specimen; multihits separate",
    original_sonic_flag=(method=="sonic"),annotation=annotation),"pooled_metrics.json",pretty=TRUE,auto_unbox=TRUE)
write_runtime("pooled_runtime.json")
# This policy controls report publication only. The original estimator is called
# with its unchanged input and minimum; rejected specimens receive no estimate.
jsonlite::write_json(list(state=report_state,abundance_failure_policy=failure_policy,
    abundance_method=method,original_sonic_minimum_bp=if(method=="sonic")
        formals(sonicLength::estAbund)$min.length else NA_integer_,
    selected_specimens=unname(unique(metadata$Specimen_ID)),
    estimated_specimens=unname(names(estimates)),unavailable_specimens=unname(failed_specimens),
    failures=unname(statuses[failed_specimens]),input_fragments_filtered=FALSE,
    unavailable_values="NA; unavailable specimens have no rows in the abundance table"),
    "original_report_status.json",pretty=TRUE,auto_unbox=TRUE)
if(length(failed_specimens)) {
    if(failure_policy=="fail")
        stop("Original estimator failed; see abundance_status.tsv. No substitute estimates were produced.")
    message("Report published with unavailable abundance for: ",paste(failed_specimens,collapse=", "),
        ". Original estimator rejection retained; no fragments excluded and no substitute estimates produced.")
}
