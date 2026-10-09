args <- commandArgs(TRUE)
suppressPackageStartupMessages(library(GenomicRanges))
read <- function(path) { e<-new.env(); vars<-load(path,e);setNames(lapply(sort(vars),function(n)e[[n]]),sort(vars)) }
checks <- list()
for(name in c("keys.RData","stats.RData","hits.R1.RData","hits.R2.RData","allSites.RData","sites.final.RData","multihitData.RData")) {
    x<-read(file.path(args[1],name));y<-read(file.path(args[2],name))
    same<-identical(x,y);checks[[name]]<-list(identical=same,detail=if(same)"" else paste(all.equal(x,y),collapse="; "))
    rm(x,y);gc()
}
x<-read(file.path(args[1],"chimeraData.RData"))$chimeraData
y<-read(file.path(args[2],"chimeraData.RData"))$chimeraData
chimera <- list(read_info_identical=identical(x$read_info,y$read_info),
    original_archived_alignments=length(x$alignments),corrected_archived_alignments=length(y$alignments),
    retained_chimera_read_rows=nrow(y$read_info),
    original_prefix_identical=identical(x$alignments,y$alignments[seq_along(x$alignments)]),
    reason="Original length(data.frame) iterates columns; nrow(data.frame) archives every retained chimera row. This does not alter insertion-site calls or chimera classification.")
chimera_ok <- chimera$read_info_identical && chimera$original_prefix_identical &&
    chimera$corrected_archived_alignments==chimera$retained_chimera_read_rows &&
    chimera$original_archived_alignments==ncol(x$read_info)
passed <- all(vapply(checks,function(z)z$identical,logical(1))) && chimera_ok
jsonlite::write_json(list(status=if(passed)"passed_with_documented_chimera_archive_fix" else "failed",
    scope="Full historical RUN809 0% prepared sequences and all 85 archived original BLAT chunks; all insertion-site, multihit, read-key, reference and caller-statistic payloads",
    checks=checks,documented_difference=chimera),args[3],pretty=TRUE,auto_unbox=TRUE)
if(!passed) stop("Historical scientific caller mismatch: ",args[3])
cat("PASS: seven complete original caller payloads identical; full chimera archive repair independently verified\n")
