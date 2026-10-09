args <- commandArgs(trailingOnly=TRUE); root <- normalizePath(args[1]); out <- args[2]
source(file.path(root,"bin/bushman_runtime.R"))
source(file.path(root,"bin/bushman_indexed.R"))
e <- bushman_runtime(file.path(root,"vendor/bushman")); fast <- bushman_optimize(e)
validate_ast <- function(x) {
    if(is.call(x)) {
        if(identical(x[[1]],as.name("<-"))) stopifnot(length(x)==3L)
        invisible(lapply(as.list(x)[-1],validate_ast))
    }
}
validate_ast(body(fast$processAlignments))
checks <- character()
for(seed in seq_len(40L)) {
    set.seed(seed)
    values <- sample(c(seq_len(400L),NA_integer_),3000L,replace=TRUE)
    queries <- IRanges::IntegerList(lapply(seq_len(100L),function(i)
        sample(c(seq_len(450L),NA_integer_),sample(0:15,1),replace=TRUE)))
    loci <- sample(seq_len(100L),500L,replace=TRUE)
    original_queries <- queries[loci]
    original <- IRanges::IntegerList(lapply(original_queries,function(q) which(values %in% q)))
    actual <- bushman_locus_membership(values,original_queries,loci)
    stopifnot(identical(actual,original))
}
stopifnot(identical(bushman_locus_membership(integer(),IRanges::IntegerList(),integer()),
    IRanges::IntegerList()))
checks <- c(checks,"randomized_locus_membership_with_NA_empty_duplicates_and_reused_loci")
keys <- data.frame(readPairKey=c("a","a","b","c","d"),ID=c("shared","one","shared","three",NA),
    medians=c(200,250,400,500,600))
components <- list(first=c("a"),second=c("b","c"),missing=c("d"))
expected <- lapply(components,function(x) {
    ids <- unique(keys[keys$readPairKey %in% x,]$ID)
    data.frame(table(keys[keys$ID %in% ids,]$medians))
})
stopifnot(identical(bushman_multihit_lengths(keys,components),expected))
for(v in list(character(),c("a"),c("a","a","b","b","b","z"))) {
    expected <- IRanges::IntegerList(lapply(unique(v),function(x) which(v==x)))
    stopifnot(identical(bushman_group_rows(v),expected))
}
checks <- c(checks,"multihit_cross_component_ID_scope_and_sorted_row_expansion")
r1 <- GenomicRanges::GRanges("chrSynthetic",IRanges::IRanges(c(100,110,150),width=10),strand="-")
r2 <- GenomicRanges::GRanges("chrSynthetic",IRanges::IRanges(c(20,30,60),width=10),strand="+")
r1$qName <- c(1L,1L,2L);r2$qName <- c(3L,3L,4L)
r1$from <- "R1";r2$from <- "R2"
reads <- data.frame(R1=c(1L,1L,2L),R2=c(3L,3L,4L),names=c("A%one","A%two","A%three"))
expected <- GenomicRanges::GRangesList(lapply(seq_len(nrow(reads)),function(i) {
    x <- r1[r1$qName==reads[i,"R1"]];y <- r2[r2$qName==reads[i,"R2"]]
    names(x) <- rep(reads[i,"names"],length(x));names(y) <- rep(reads[i,"names"],length(y))
    c(y,x)
}))
stopifnot(identical(bushman_chimera_alignments(r1,r2,reads),expected))
stopifnot(identical(bushman_chimera_alignments(r1,r2,reads[FALSE,]),GenomicRanges::GRangesList()))
checks <- c(checks,"chimera_all_hits_original_multiplicity_names_and_R2_R1_order")
dir.create(file.path(out,"parser"),recursive=TRUE,showWarnings=FALSE)
write_gzip <- function(path,lines) { con <- gzfile(path,"wt");writeLines(lines,con);close(con);path }
line <- paste(c(97,3,0,0,0,0,0,5,"+","001",100,5,100,"chrSynthetic",24000,1000,1095,1,"95,","5,","1000,"),collapse="\t")
paths <- c(write_gzip(file.path(out,"parser/empty.psl.gz"),character()),
    write_gzip(file.path(out,"parser/one.psl.gz"),line),
    write_gzip(file.path(out,"parser/two.psl.gz"),c(line,sub("001","02",line,fixed=TRUE))))
original_env <- bushman_runtime(file.path(root,"vendor/bushman"))
for(p in list(character(),paths[1],paths[2],paths[c(3,1,2)])) {
    stopifnot(identical(bushman_readpsl_stream(p),original_env$readpsl(p)))
    stopifnot(identical(bushman_readpsl_stream(p,c("qStarts","tStarts")),original_env$readpsl(p,c("qStarts","tStarts"))))
}
checks <- c(checks,"streaming_PSL_types_empty_single_and_ordered_multi_file_parity","five_frozen_caller_anchors")
jsonlite::write_json(list(status="passed",checks=checks,random_seeds=40L),file.path(out,"kernel_results.json"),
    pretty=TRUE,auto_unbox=TRUE)
cat("PASS: indexed kernels and streaming parser preserve baseline payloads\n")
