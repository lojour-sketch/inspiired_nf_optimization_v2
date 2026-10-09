args <- commandArgs(TRUE)
source(file.path(args[1],"vendor/bushman/intSiteCaller/linker_common.R"))
values <- jsonlite::fromJSON(args[2])
original <- linker_common(values$linkerSequence)
stopifnot(identical(unname(original),values$python_effective))
jsonlite::write_json(list(status="passed",oracle="committed upstream linker_common.R",
    cases=values,original_effective=unname(original),
    RUN809_effective_marker="AGTCCCTTAAGCGGAG",raw_CSV_forward_marker="CTCCGCTTAAGGGACT"),
    args[3],pretty=TRUE,auto_unbox=TRUE)
cat("PASS: Python input normalization matches original effective metadata\n")
