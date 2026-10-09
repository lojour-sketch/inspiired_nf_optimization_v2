args <- commandArgs(TRUE)
suppressPackageStartupMessages(library(data.table))
original <- fread(args[1],data.table=FALSE)
corrected <- fread(args[2],data.table=FALSE)
fields <- c("seqnames","start","strand","posid","estAbund","estAbundProp","estAbundRank")
original <- original[order(original$posid),fields]
corrected <- corrected[order(corrected$posid),fields]
rownames(original)<-rownames(corrected)<-NULL
# Text round-trips can introduce insignificant floating-point differences.
exact <- identical(original[,setdiff(fields,"estAbundProp")],corrected[,setdiff(fields,"estAbundProp")])
proportion <- isTRUE(all.equal(original$estAbundProp,corrected$estAbundProp,tolerance=1e-12))
passed <- exact && proportion
jsonlite::write_json(list(status=if(passed)"passed" else "failed",
    scope="Historical original 0% PCR evidence recomputed with frozen original 5 bp/Sonic reporting; corrected single-specimen report at the same scope",
    coordinate_abundance_rank_exact=exact,proportions_equal_tolerance_1e_12=proportion,
    original_sites=nrow(original),corrected_sites=nrow(corrected),
    original_rounded_abundance=sum(original$estAbund),corrected_rounded_abundance=sum(corrected$estAbund),
    largest_unique_site_fraction=max(corrected$estAbundProp),
    details=if(passed)"" else paste(all.equal(original,corrected),collapse="; ")),
    args[3],pretty=TRUE,auto_unbox=TRUE)
if(!passed) stop("Original single-specimen Sonic report mismatch: ",args[3])
cat("PASS: original coordinates, rounded Sonic abundance, proportions and ranks recovered\n")
