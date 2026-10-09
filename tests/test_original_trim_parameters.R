args <- commandArgs(trailingOnly=TRUE)
pipeline <- normalizePath(args[1]);out <- normalizePath(args[2])
source(file.path(pipeline,"bin/bushman_runtime.R"))
e <- bushman_runtime(file.path(pipeline,"vendor/bushman"))
marker <- "CTCCGCTTAAGGGACT";prefix <- "ACGTCGATCGTTAGC"
reads <- DNAStringSet(c(with_marker=paste0(prefix,marker,"TAACCG"),without_marker=strrep("A",50)))
expected <- c(nchar(prefix),50-nchar(marker)/2)
result <- e$trim_overreading(reads,marker,3)
stopifnot(identical(as.numeric(width(result)),as.numeric(expected)))
reversed <- e$trim_overreading(reads,as.character(reverseComplement(DNAString(marker))),3)
stopifnot(!identical(width(result),width(reversed)))
jsonlite::write_json(list(status="passed",marker=marker,rule="Supplied historical orientation and original batch-dependent precautionary clipping",
    input_widths=as.numeric(width(reads)),expected_widths=expected,observed_widths=as.numeric(width(result)),
    reverse_complement_changes_result=TRUE),file.path(out,"original_trim_parameter_check.json"),pretty=TRUE,auto_unbox=TRUE)
cat("PASS: historical read-through orientation and batch rule\n")
