process LTRchecking {

    memory '80 GB'

    input:
    tuple val(meta), path(read1), path(read2), val(primer), val(ltrbit), val(largeLTRFrag), val(idx)

    output:
    tuple val(meta.id), path("${meta.id}.ltr_filtered_R1.fastq.gz"), path("${meta.id}.ltr_filtered_R2.fastq.gz"), path("${meta.id}.ltr_checked.log")


    script:
    """
    #!/usr/bin/env Rscript

    library(Biostrings)
    library(ShortRead)

    source("/home/lrenteria/inspiired_nf/usefulscripts/pairwiseAlignment.R")

    primer <- "${primer}"
    ltrbit <- "${ltrbit}"

    ### substitution matrix
    submat <- nucleotideSubstitutionMatrix(match = 1, mismatch = 0, baseOnly = TRUE)

    #### assign variables and read files
    r1_file <- "${read1}"
    r2_file <- "${read2}"
    out_r1 <- "${meta.id}.ltr_filtered_R1.fastq.gz"
    out_r2 <- "${meta.id}.ltr_filtered_R2.fastq.gz"
    log_file <- "${meta.id}.ltr_checked.log"

    reads_r1 <- readFastq(r1_file)
    reads_r2 <- readFastq(r2_file)

    seqs_r2 <- sread(reads_r2)

    ### primer alignment
    aln_p <- pairwiseAlignment(
        pattern = subseq(seqs_r2, 1, pmin(width(seqs_r2), 1 + nchar(primer))),
        subject = primer,
        substitutionMatrix = submat,
        gapOpening = 0,
        gapExtension = 1,
        type = "overlap"
    )
    aln_p_df <- PairwiseAlignmentsSingleSubject2DF(aln_p)

    ### ltr alignment
    start_ltr <- nchar(primer) + 1
    end_ltr <- start_ltr + nchar(ltrbit) + 5
    aln_l <- pairwiseAlignment (
        pattern = subseq(seqs_r2, start_ltr, pmin(width(seqs_r2), end_ltr)),
        subject = ltrbit,
        substitutionMatrix = submat,
        gapOpening = 0,
        gapExtension = 1,
        type = "overlap"
    )
    aln_l_df <- PairwiseAlignmentsSingleSubject2DF(aln_l, shift = nchar(primer))

    ### filtering by score
    max_mismatches <- 3
    primer_ok <- aln_p_df\$nmismatch <= max_mismatches
    ltr_ok    <- aln_l_df\$nmismatch <= max_mismatches
    keep_mask <- primer_ok & ltr_ok
    keep_ids <- as.character(id(reads_r2))[keep_mask]

    ### create output
    writeFastq(reads_r1[id(reads_r1) %in% keep_ids], out_r1, compress = TRUE)
    writeFastq(reads_r2[id(reads_r2) %in% keep_ids], out_r2, compress = TRUE)

    cat("Sample: ${meta.id}\n",
    "Total R2 reads: ", length(reads_r2), "\n",
    "Reads kept: ", length(keep_ids), "\n",
    file = log_file)
    """


}