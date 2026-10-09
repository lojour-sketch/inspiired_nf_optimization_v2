args <- commandArgs(trailingOnly=TRUE)
if(length(args)!=5L) stop("Usage: BASELINE_RESULTS BASELINE_PROJECT FAST_RESULTS FAST_PROJECT OUTPUT_JSON")
suppressPackageStartupMessages({library(GenomicRanges);library(data.table)})
base <- args[1];fast <- args[3];out <- args[5]
results <- list()
compare <- function(label,x,y) {
    equal <- identical(x,y)
    detail <- if(equal) "" else paste(all.equal(x,y,check.attributes=TRUE),collapse="; ")
    results[[length(results)+1L]] <<- list(item=label,identical=equal,detail=detail)
}
load_payload <- function(path) {
    env <- new.env(parent=emptyenv()); names <- load(path,envir=env)
    setNames(lapply(sort(names),function(n) {
        x <- env[[n]]
        # XStringSet's shared storage contains external pointers. Compare its
        # complete scientific contents, not serialized memory addresses.
        if(methods::is(x,"XStringSet")) return(list(class=class(x),sequences=as.character(x),
            widths=Biostrings::width(x),names=names(x),metadata=S4Vectors::metadata(x),
            mcols=as.data.frame(S4Vectors::mcols(x))))
        x
    }),sort(names))
}
calls1 <- file.path(base,"04_original_calls",args[2]);calls2 <- file.path(fast,"04_original_calls",args[4])
samples1 <- sort(list.dirs(calls1,recursive=FALSE,full.names=FALSE))
samples2 <- sort(list.dirs(calls2,recursive=FALSE,full.names=FALSE))
if(!length(samples1) || !identical(samples1,samples2)) stop("Missing or unexpected samples")
files <- c("keys.RData","stats.RData","primerIDData.RData","hits.R1.RData","hits.R2.RData","allSites.RData",
    "sites.final.RData","multihitData.RData","chimeraData.RData")
for(sid in samples1) for(name in files) {
    a <- file.path(calls1,sid,sid,name);b <- file.path(calls2,sid,sid,name)
    if(!file.exists(a) && !file.exists(b) && name %in% c("primerIDData.RData","hits.R1.RData","hits.R2.RData")) next
    if(!file.exists(a) || !file.exists(b)) stop("Missing payload: ",sid,"/",name)
    compare(paste(sid,name,sep="/"),load_payload(a),load_payload(b))
}
for(sid in samples1) {
    a <- file.path(calls1,sid,sid);b <- file.path(calls2,sid,sid)
    fa <- sort(list.files(a,pattern="^R[12]-.*\\.fa$"))
    if(!identical(fa,sort(list.files(b,pattern="^R[12]-.*\\.fa$")))) stop("Curated FASTA chunk mismatch: ",sid)
    for(name in fa) compare(paste(sid,name,sep="/"),readLines(file.path(a,name)),readLines(file.path(b,name)))
}
pool1 <- file.path(base,"05_original_report",args[2]);pool2 <- file.path(fast,"05_original_report",args[4])
groups1 <- sort(list.dirs(pool1,recursive=FALSE,full.names=FALSE))
groups2 <- sort(list.dirs(pool2,recursive=FALSE,full.names=FALSE))
if(!length(groups1) || !identical(groups1,groups2)) stop("Missing or unexpected pooled report groups")
for(group in groups1) {
    a <- file.path(pool1,group);b <- file.path(pool2,group)
    payloads <- sort(list.files(a,pattern="\\.rds$|\\.tsv$"))
    for(name in payloads) {
        if(!file.exists(file.path(b,name))) stop("Missing pooled payload: ",name)
        x <- if(grepl("\\.rds$",name)) readRDS(file.path(a,name)) else fread(file.path(a,name),data.table=FALSE)
        y <- if(grepl("\\.rds$",name)) readRDS(file.path(b,name)) else fread(file.path(b,name),data.table=FALSE)
        compare(paste(group,name,sep="/"),x,y)
    }
}
passed <- all(vapply(results,function(x) x$identical,logical(1)))
jsonlite::write_json(list(status=if(passed) "passed" else "failed",baseline=base,optimized=fast,
    samples=samples1,checks=results),out,pretty=TRUE,auto_unbox=TRUE)
if(!passed) stop("Scientific payload mismatch; inspect ",out)
cat("PASS: ",length(results)," complete caller and pooled payload comparisons\n",sep="")
