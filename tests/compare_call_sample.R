args <- commandArgs(trailingOnly=TRUE)
suppressPackageStartupMessages(library(GenomicRanges))
results <- list()
for(name in c("keys.RData","stats.RData","hits.R1.RData","hits.R2.RData",
    "allSites.RData","sites.final.RData","multihitData.RData","chimeraData.RData")) {
    read <- function(path) {
        env <- new.env(parent=emptyenv());vars <- load(path,env)
        setNames(lapply(sort(vars),function(n) env[[n]]),sort(vars))
    }
    x <- read(file.path(args[1],name));y <- read(file.path(args[2],name))
    equal <- identical(x,y)
    results[[name]] <- list(identical=equal,detail=if(equal) "" else paste(all.equal(x,y),collapse="; "))
    rm(x,y);invisible(gc())
}
passed <- all(vapply(results,function(x) x$identical,logical(1)))
jsonlite::write_json(list(status=if(passed) "passed" else "failed",checks=results),args[3],
    pretty=TRUE,auto_unbox=TRUE)
if(!passed) stop("Real-data caller mismatch: ",args[3])
cat("PASS: all eight real-data caller archives exactly match the baseline\n")
