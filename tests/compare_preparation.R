args <- commandArgs(TRUE)
suppressPackageStartupMessages(library(GenomicRanges))
checks <- list()
compare <- function(name,x,y) {
    equal <- identical(x,y)
    checks[[name]] <<- list(identical=equal,detail=if(equal)"" else paste(all.equal(x,y),collapse="; "))
}
read <- function(path) { e <- new.env();n <- load(path,e);e[[n[1]]] }
columns <- c("sample","barcoded","LTRed","linkered","ltredlinkered","lenTrim","vTrimed","uniqL","uniqP","uniqP30","lLen","pLen")
x <- read(file.path(args[1],"stats.RData"));y <- read(file.path(args[2],"stats.RData"))
compare("trim_statistics",x[,columns],y[,columns])
statistics <- list(original=x[,columns],corrected=y[,columns])
rm(x,y);gc()
for(name in c("keys.RData","primerIDData.RData")) {
    x <- read(file.path(args[1],name));y <- read(file.path(args[2],name))
    canonical <- function(z) {
        if(methods::is(z,"XStringSet")) return(list(class=class(z),sequence=as.character(z),
            names=names(z),mcols=as.data.frame(mcols(z)),metadata=metadata(z)))
        z
    }
    compare(name,canonical(x),canonical(y));rm(x,y);gc()
}
files <- sort(list.files(args[1],pattern="^R[12]-.*\\.fa$"))
compare("FASTA_chunk_names",files,sort(list.files(args[2],pattern="^R[12]-.*\\.fa$")))
for(name in files) {
    compare(name,readLines(file.path(args[1],name)),readLines(file.path(args[2],name)))
    gc()
}
passed <- all(vapply(checks,function(x)x$identical,logical(1)))
jsonlite::write_json(list(status=if(passed)"passed" else "failed",checks=checks,statistics=statistics,
    scope="Full RUN809 0%: all curated sequence text/chunks, paired-read keys and complete primer-ID contents; upstream effective linker marker"),
    args[3],pretty=TRUE,auto_unbox=TRUE)
if(!passed) stop("Preparation differs from historical original: ",args[3])
cat("PASS: complete original preparation recovered after effective linker correction\n")
