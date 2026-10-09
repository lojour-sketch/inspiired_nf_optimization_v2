args <- commandArgs(trailingOnly=TRUE)
root <- normalizePath(args[1]); out <- normalizePath(args[2])
source(file.path(root,"bin/bushman_runtime.R"))
source(file.path(root,"bin/bushman_indexed.R"))
source(file.path(root,"bin/bushman_compact.R"))
native <- bushman_load_pairing(file.path(root,"bin/bushman_pairing.c"))
e <- bushman_compact_optimize(bushman_runtime(file.path(root,"vendor/bushman")),native)
stopifnot(length(e$optimization_contract$compact_anchors)==7L)
# Candidate geometry keeps names, ranges, strand and complete sequence dictionary,
# while the large revmap remains only on the reduced loci used for membership.
toy <- GenomicRanges::GRanges(c("chrA","chrB"),IRanges::IRanges(c(10,50),width=c(2,3)),strand=c("+","-"))
toy$revmap <- IRanges::IntegerList(list(1:1000,1001:2000))
ids <- c(2L,1L,2L,2L)
stopifnot(identical(GenomicRanges::granges(toy[ids]),GenomicRanges::granges(toy,use.mcols=FALSE)[ids]))
oracle <- function(keys,h1,h2,r1,r2,loci) {
    queries <- function(h,r,ids) lapply(r$revmap[ids],function(x) as.integer(h$qName[x]))
    a <- lapply(queries(h1,r1,loci$R1.loci),function(q) which(keys$R1 %in% q))
    b <- lapply(queries(h2,r2,loci$R2.loci),function(q) which(keys$R2 %in% q))
    IRanges::CharacterList(lapply(seq_len(nrow(loci)),function(i)
        as.character(keys$readPairKey[intersect(a[[i]],b[[i]])])))
}
checks <- character()
for(seed in seq_len(60L)) {
    set.seed(seed)
    key_values <- c(as.character(1:80),"001",NA_character_)
    keys <- data.frame(R1=sample(key_values,1500,replace=TRUE),
        R2=sample(key_values,1500,replace=TRUE),readPairKey=sample(paste0("pair",1:900),1500,replace=TRUE))
    h1 <- list(qName=sample(c(as.character(1:95),"001",NA_character_),500,replace=TRUE))
    h2 <- list(qName=sample(c(as.character(1:95),"002",NA_character_),500,replace=TRUE))
    r1 <- list(revmap=IRanges::IntegerList(lapply(1:50,function(i) sample(1:500,sample(0:30,1),replace=TRUE))))
    r2 <- list(revmap=IRanges::IntegerList(lapply(1:40,function(i) sample(1:500,sample(0:30,1),replace=TRUE))))
    loci <- data.frame(R1.loci=sample(1:50,400,replace=TRUE),R2.loci=sample(1:40,400,replace=TRUE))
    expected <- oracle(keys,h1,h2,r1,r2,loci)
    actual <- bushman_compact_support(keys,h1,h2,r1,r2,loci,native,function(x) {})
    stopifnot(identical(expected,actual))
    # Small blocks exercise a different join partition without changing support.
    a <- bushman_compact_members(keys$R1,h1,r1,loci$R1.loci,block_size=3L)
    expected_members <- IRanges::IntegerList(lapply(r1$revmap[loci$R1.loci],
        function(x) which(keys$R1 %in% as.integer(h1$qName[x]))))
    got <- relist(a$values,IRanges::PartitioningByEnd(a$ends))[a$slots]
    stopifnot(identical(got,expected_members))
    empty <- loci[FALSE,]
    stopifnot(identical(bushman_compact_support(keys,h1,h2,r1,r2,empty,native,function(x) {}),
        IRanges::CharacterList()))
}
checks <- c(checks,"60 randomized original scan/intersection oracles: duplicate rows/IDs/queries, NA, leading zero queries, missing keys, shuffled/reused loci, empty inputs and block boundaries")
for(seed in 1:100) {
    set.seed(seed)
    a <- sort(sample(1:5000,sample(0:200,1)))
    b <- sort(sample(1:5000,sample(0:4500,1)))
    got <- .Call(native$symbol,a,as.integer(length(a)),b,as.integer(length(b)),1L,1L)
    stopifnot(identical(got$values,intersect(a,b)),identical(got$ends,as.integer(length(intersect(a,b)))))
}
fails <- function(...) inherits(try(.Call(native$symbol,...),silent=TRUE),"try-error")
stopifnot(fails(c(2L,1L),2L,1L,1L,1L,1L),fails(c(1L,1L),2L,1L,1L,1L,1L),
    fails(NA_integer_,1L,1L,1L,1L,1L),fails(1L,2L,1L,1L,1L,1L),
    fails(1L,1L,1L,1L,2L,1L),fails(1L,1L,1L,1L,NA_integer_,1L))
checks <- c(checks,"100 randomized balanced/unequal native intersections and invalid-input guards",
    "seven compact AST anchors plus five existing indexed expansion anchors; candidate geometry unchanged without repeated revmap")

# Dense hub microbenchmark: same stored input and identical complete support;
# original and old indexed implementations both repeat the large R2 membership.
n <- 80000L; nl <- 240L
keys <- data.frame(R1=rep(seq_len(nl),length.out=n),R2=rep(1L,n),readPairKey=paste0("k",seq_len(n)))
h1 <- list(qName=seq_len(nl));h2 <- list(qName=1L)
r1 <- list(revmap=IRanges::IntegerList(as.list(seq_len(nl))))
r2 <- list(revmap=IRanges::IntegerList(list(1L)))
loci <- data.frame(R1.loci=seq_len(nl),R2.loci=rep(1L,nl))
run_original <- function() oracle(keys,h1,h2,r1,r2,loci)
run_indexed <- function() {
    q1 <- IRanges::IntegerList(as.list(seq_len(nl)))
    q2 <- IRanges::IntegerList(rep(list(1L),nl))
    a <- bushman_locus_membership(keys$R1,q1,loci$R1.loci)
    b <- bushman_locus_membership(keys$R2,q2,loci$R2.loci)
    IRanges::CharacterList(lapply(seq_len(nl),function(i) keys$readPairKey[intersect(a[[i]],b[[i]])]))
}
run_compact <- function() bushman_compact_support(keys,h1,h2,r1,r2,loci,native,function(x) {})
expected <- run_original()
stopifnot(identical(expected,run_indexed()),identical(expected,run_compact()))
timings <- list()
for(name in c("original","indexed","compact")) {
    fn <- get(paste0("run_",name))
    timings[[name]] <- vapply(1:3,function(i) {gc();unname(system.time(result <- fn())["elapsed"])},numeric(1))
}
summary <- list(status="passed",checks=checks,random_seeds=60L,native=native$provenance,
    benchmark=list(scope="dense-hub pairing kernel only; excludes compilation, BLAT, caller classifications and report",
        key_rows=n,candidate_pairs=nl,distinct_R1_loci=nl,distinct_R2_loci=1L,
        elapsed_seconds=timings,median_seconds=lapply(timings,median),
        indexed_membership_rows=n*nl+n,compact_membership_rows=2L*n,
        support_rows=n,complete_support_identical=TRUE))
jsonlite::write_json(summary,file.path(out,"compact_kernel_results.json"),pretty=TRUE,auto_unbox=TRUE)
cat("PASS: compact memberships/intersections match original operations; dense-hub benchmark recorded\n")
