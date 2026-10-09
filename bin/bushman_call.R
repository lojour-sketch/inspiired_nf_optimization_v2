#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly=TRUE)
options(error=function() { traceback(3);quit(status=1) })
if(!length(args) %in% c(5L,6L)) stop("Usage: ROW_JSON SAMPLE_DIRECTORY TWOBIT SOURCES PSL_DIRECTORY [compact|indexed|original]")
own <- dirname(normalizePath(sub("^--file=","",grep("^--file=",commandArgs(),value=TRUE)[1])))
task_dir <- getwd()
source(file.path(own,"bushman_runtime.R"))
row <- jsonlite::fromJSON(args[1])
e <- bushman_runtime(args[4],reference=args[3],reference_seqinfo=row$reference_seqinfo_rds)
implementation <- if(length(args)==6L) args[6] else "indexed"
if(!implementation %in% c("compact","indexed","original")) stop("Unknown caller implementation")
if(implementation %in% c("compact","indexed")) {
    source(file.path(own,"bushman_indexed.R"))
    if(implementation=="compact") {
        source(file.path(own,"bushman_compact.R"))
        native <- bushman_load_pairing(file.path(own,"bushman_pairing.c"))
        e <- bushman_compact_optimize(e,native)
    } else e <- bushman_optimize(e)
}
original_mark <- e$mark
e$mark <- function(label) {
    original_mark(label)
    jsonlite::write_json(list(sample=row$Sample_ID,implementation=implementation,
        phase=label,updated_at_utc=format(Sys.time(),tz="UTC",usetz=TRUE),
        completed_phase_seconds=e$timings),file.path(task_dir,"calling_progress.json"),
        pretty=TRUE,auto_unbox=TRUE)
}
sample <- row$Sample_ID
if(!dir.exists(sample)) dir.create(sample)
file.copy(list.files(args[2],full.names=TRUE),sample,overwrite=TRUE,recursive=TRUE)
psl <- list.files(args[5],pattern="\\.psl\\.gz$",full.names=TRUE)
file.copy(psl,sample,overwrite=TRUE)
codeDir <- normalizePath(file.path(args[4],"intSiteCaller"))
save(codeDir,file="codeDir.RData")
keys <- get(load(file.path(sample,"keys.RData")))
started <- proc.time()
profile <- identical(Sys.getenv("BUSHMAN_PROFILE"),"1")
if(profile) Rprof(file.path(task_dir,"calling.Rprof"),memory.profiling=TRUE,gc.profiling=TRUE)
if(nrow(keys)) {
    expected <- basename(list.files(sample,pattern="^R[12]-.*\\.fa$"))
    if(!setequal(paste0(expected,".psl.gz"),basename(psl))) stop("Incomplete BLAT chunk set")
    e$processAlignments(sample,as.numeric(row$minPctIdent),as.integer(row$maxAlignStart),
        as.integer(row$maxFragLength),row$refGenome)
} else {
    allSites <- sites.final <- GenomicRanges::GRanges()
    multihitData <- list(unclusteredMultihits=GenomicRanges::GRanges(),
        clusteredMultihitPositions=GenomicRanges::GRangesList(),clusteredMultihitLengths=list())
    chimeraData <- list(read_info=keys,alignments=GenomicRanges::GRangesList())
    save(allSites,file=file.path(sample,"allSites.RData"))
    save(sites.final,file=file.path(sample,"sites.final.RData"))
    save(multihitData,file=file.path(sample,"multihitData.RData"))
    save(chimeraData,file=file.path(sample,"chimeraData.RData"))
}
setwd(task_dir)
if(profile) Rprof(NULL)
e$mark("finished")
jsonlite::write_json(list(sample=sample,input_curated_pairs=nrow(keys),psl_chunks=length(psl),
    elapsed_seconds=unname((proc.time()-started)[3]),
    phase_seconds=e$timings,
    kernel_seconds=e$kernel_seconds,compact_pairing_metrics=e$compact_pairing_metrics,
    implementation=implementation,
    optimization_contract=e$optimization_contract,
    algorithm="original independent BLAT mates, all qualifying locus combinations, original filters and multihit graph"),
    "calling_metrics.json",pretty=TRUE,auto_unbox=TRUE)
write_runtime("calling_runtime.json")
runtime <- jsonlite::fromJSON("calling_runtime.json",simplifyVector=FALSE)
identity <- jsonlite::fromJSON(file.path(own,"..","pipeline_identity.json"))
runtime$pipeline <- identity$pipeline;runtime$version <- identity$version
runtime$caller_implementation <- implementation
runtime$optimization <- e$optimization_contract
jsonlite::write_json(runtime,"calling_runtime.json",pretty=TRUE,auto_unbox=TRUE)
