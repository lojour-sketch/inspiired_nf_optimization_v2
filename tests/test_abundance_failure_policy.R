#!/usr/bin/env Rscript
args <- commandArgs(TRUE)
exprs <- parse(args[1])
last <- exprs[[length(exprs)]]
# Exercise the actual publication guard with a real rejection from Sonic 1.4.7.
rejection <- tryCatch(sonicLength::estAbund(c("site1","site1"),c(2L,30L)),
    error=function(err) conditionMessage(err))
stopifnot(is.character(rejection),grepl("min.length",rejection,fixed=TRUE),
    identical(as.numeric(formals(sonicLength::estAbund)$min.length),20))
run <- function(policy,failed) {
    env <- new.env(parent=globalenv());env$failure_policy <- policy;env$failed_specimens <- failed
    tryCatch({eval(last,env);"published"},error=function(err) conditionMessage(err))
}
stopifnot(grepl("Original estimator failed",run("fail","Polyclonal_95pct"),fixed=TRUE),
    identical(run("record","Polyclonal_95pct"),"published"),
    identical(run("fail",character()),"published"),identical(run("record",character()),"published"))
# Seven-argument callers retain the strict default.
initial <- new.env();initial$args <- rep("unused",7L)
eval(exprs[[3]],initial);stopifnot(identical(initial$failure_policy,"fail"))
initial$args <- c(rep("unused",7L),"record");eval(exprs[[3]],initial)
stopifnot(identical(initial$failure_policy,"record"))
jsonlite::write_json(list(status="passed",Sonic_minimum_bp=20L,
    original_short_fragment_rejection=rejection,default_policy="fail",explicit_record_policy="publishes missing abundance without an estimate",
    cases=c("strict rejection","explicit partial publication","valid strict publication","valid record publication","legacy strict default")),
    args[2],pretty=TRUE,auto_unbox=TRUE)
cat("PASS: original Sonic minimum and rejection retained; explicit publication policy verified\n")
