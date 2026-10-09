#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly=TRUE)
if(length(args)!=7L) stop("Usage: ROW_JSON R1 R2 VECTOR SOURCES BLAT CHUNK_SIZE")
own <- dirname(normalizePath(sub("^--file=","",grep("^--file=",commandArgs(),value=TRUE)[1])))
task_dir <- getwd()
r1_file <- normalizePath(args[2],mustWork=TRUE)
r2_file <- normalizePath(args[3],mustWork=TRUE)
vector_file <- normalizePath(args[4],mustWork=TRUE)
vector_name <- "bushman_vector_input.fa"
file.copy(vector_file,vector_name,overwrite=TRUE)
source(file.path(own,"bushman_runtime.R"))
e <- bushman_runtime(args[5]); row <- jsonlite::fromJSON(args[1])
Sys.setenv(PATH=paste(dirname(normalizePath(args[6])),Sys.getenv("PATH"),sep=":"))
# Synchronization is an input contract, never an alignment/abundance filter.
r1 <- ShortRead::readFastq(args[2]); r2 <- ShortRead::readFastq(args[3])
ids <- function(x) sub(" .*$","",as.character(ShortRead::id(x)))
r1_ids <- ids(r1);r2_ids <- ids(r2)
if(length(r1)!=length(r2) || !identical(r1_ids,r2_ids) || anyDuplicated(r1_ids))
    stop("FASTQ pairs must have synchronized, unique read IDs")
n_input <- length(r1)
# Reuse the very same validated objects in the upstream loader. Its whole-
# sample trimming decisions, names, sequences and quality data are unchanged.
e$validated_reads <- list(setNames(list(r1),r1_file),setNames(list(r2),r2_file))
e$consume_validated_reads <- function(read1,read2) {
    if(!identical(read1,r1_file) || !identical(read2,r2_file)) stop("Validated FASTQ paths changed")
    result <- e$validated_reads;e$validated_reads <- NULL;result
}
statements <- as.list(body(e$getTrimmedSeqs));replaced <- 0L
for(i in seq_along(statements)) {
    statement <- statements[[i]]
    if(is.call(statement) && identical(statement[[1]],as.name("<-")) &&
       identical(statement[[2]],as.name("reads")) && is.call(statement[[3]]) &&
       identical(statement[[3]][[1]],as.name("lapply"))) {
        statement[[3]] <- quote(consume_validated_reads(read1,read2))
        statements[[i]] <- statement;replaced <- replaced+1L
    }
}
if(replaced!=1L) stop("Frozen FASTQ loader anchor changed")
body(e$getTrimmedSeqs) <- as.call(statements)
rm(r1,r2,r1_ids,r2_ids);invisible(gc())
codeDir <- normalizePath(file.path(args[5],"intSiteCaller"))
save(codeDir,file="codeDir.RData")
e$config <- list(chunkSize=as.integer(args[7]))
dir.create(row$Sample_ID)
started <- proc.time()
if(n_input) {
    outcome <- tryCatch(e$getTrimmedSeqs(qualityThreshold=as.numeric(row$qualityThreshold),
        badQuality=as.numeric(row$badQualityBases),qualityWindow=as.numeric(row$qualitySlidingWindow),
        primer=row$primer,ltrbit=row$ltrbit,largeLTRFrag=row$largeLTRFrag,
        linker=row$linkerSequence,linker_common=row$linkerCommon,mingDNA=as.integer(row$mingDNA),
        read1=r1_file,read2=r2_file,alias=row$Sample_ID,
        vectorSeq=vector_name),error=function(err) {
            if(identical(conditionMessage(err),"error - no curated reads")) return("empty")
            stop(err)
        })
} else outcome <- "empty"
setwd(task_dir)
e$mark("finished")
sample_dir <- row$Sample_ID
if(!file.exists(file.path(sample_dir,"keys.RData"))) {
    keys <- data.frame(R2=integer(),R1=integer(),names=character(),readPairKey=character())
    save(keys,file=file.path(sample_dir,"keys.RData"))
    stats <- data.frame(sample=row$Sample_ID,barcoded=n_input)
    save(stats,file=file.path(sample_dir,"stats.RData"))
}
file.copy("codeDir.RData",file.path(sample_dir,"codeDir.RData"),overwrite=TRUE)
jsonlite::write_json(list(sample=row$Sample_ID,input_pairs=n_input,
    accepted_pairs=nrow(get(load(file.path(sample_dir,"keys.RData")))),
    elapsed_seconds=unname((proc.time()-started)[3]),chunk_size=as.integer(args[7]),
    phase_seconds=e$timings,
    logic="frozen Bushman trimming/vector filtering; no Trim Galore, regex linker filter or minimap2"),
    "trim_metrics.json",pretty=TRUE,auto_unbox=TRUE)
write_runtime("trim_runtime.json")
chunks <- list.files(sample_dir,pattern="^R[12]-.*\\.fa$",full.names=TRUE)
expected_psls <- if(length(chunks)) paste0(basename(chunks),".psl.gz") else character()
jsonlite::write_json(list(sample=row$Sample_ID,expected_psl_names=expected_psls),
    file.path(sample_dir,"chunks_manifest.json"),pretty=TRUE,auto_unbox=FALSE)
write.table(data.frame(sample=rep(row$Sample_ID,length(chunks)),file=chunks),"chunks.tsv",sep="\t",quote=FALSE,row.names=FALSE)
