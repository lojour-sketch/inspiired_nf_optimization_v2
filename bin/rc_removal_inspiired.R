#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)

# Parse arguments
meta <- args[1]
read1 <- args[2] 
read2 <- args[3]
largeLTRfrag <- args[4]
common_linker <- args[5]

# DEBUGGING - Print everything
cat("=== DEBUG INFO ===\n")
cat("meta:", meta, "\n")
cat("read1:", read1, "\n")
cat("read2:", read2, "\n")
cat("read1 exists:", file.exists(read1), "\n")
cat("read2 exists:", file.exists(read2), "\n")
cat("read1 size:", file.info(read1)$size, "bytes\n")
cat("read2 size:", file.info(read2)$size, "bytes\n")
cat("Working directory:", getwd(), "\n")
cat("Files in working directory:\n")
print(list.files())
cat("==================\n")

library(Biostrings)
library(ShortRead)
library(pwalign) #we need it for the pairwiseAlignment and nucleotideSubstitutionMatrix functions

# Resolve to absolute real paths
read1 <- normalizePath(read1, mustWork = TRUE)
read2 <- normalizePath(read2, mustWork = TRUE)

fq1 <- readFastq(read1)
fq2 <- readFastq(read2)

if (length(fq1) == 0 || length(fq2) == 0) {
    cat("WARNING: One or both input files are empty. Creating empty output files. This could be due to the processing of control samples, or due to a problem with the input files.\n")
    
    # Create empty fastq files
    output_file1 <- paste0(meta, ".rc_removed_R1.fastq.gz")
    output_file2 <- paste0(meta, ".rc_removed_R2.fastq.gz")
    
    # Write empty compressed fastq files
    empty_fq <- ShortReadQ()
    writeFastq(empty_fq, output_file1, compress=TRUE)
    writeFastq(empty_fq, output_file2, compress=TRUE)
    
    cat("Created empty output files successfully.\n")
    quit(status=0)  # Exit successfully
}

# Verify they have the same number of reads initially
if(length(fq1) != length(fq2)) {
  stop("Initial read count mismatch: R1=", length(fq1), " R2=", length(fq2))
}

#R needs to convert the quality and the sequence to an specific format
reads1 <- sread(fq1)
reads2 <- sread(fq2)

if (length(reads1) == 0) stop("reads1 has length 0")
if (length(reads2) == 0) stop("reads2has length 0")

qual1 <- quality(fq1)
qual2 <- quality(fq2)

#we obtain the reverse complements of the markers
largeLTRfrag_rc <- reverseComplement(DNAString(largeLTRfrag))
common_linker_rc <- reverseComplement(DNAString(common_linker))

# we first make a function that converts the output of the pairwise alignment into a dataframe
PairwiseAlignmentsSingleSubject2DF <- function(PA, shift=0) {
    return(data.frame(
        width=width(pattern(PA)),
        score=score(PA),
        mismatch=width(pattern(PA))-score(PA),
        start=start(pattern(PA))+shift,
        end=end(pattern(PA))+shift
    ))
}

#now our function to perform the alignment and apply the INSPIIRED criteria
#we will also add the quality variable to the inputs, because we need the quality string to be cut exactly as the sequence
#inspiired does not do this because it has a specific cell in its metadata table for the sequence data, and it doesnt depend on fastq writing
trim_overreading <- function(reads, qual, marker, maxMismatch=3) {

    submat <- pwalign::nucleotideSubstitutionMatrix(match=1, mismatch=0, baseOnly=TRUE)

    #allows gap, and del/ins count as 1 mismatch
    tmp <- pwalign::pairwiseAlignment(pattern=reads,
                                subject=marker,
                                substitutionMatrix=submat,
                                gapOpening=0,
                                gapExtension=1,
                                type="overlap")
    
    if (length(pattern(tmp)) == 0) stop("pairwiseAlignment returned zero patterns")

    #we convert the pairwisealignment output into a dataframe
    odf <- PairwiseAlignmentsSingleSubject2DF(tmp)

    #we assume no reads have good alignments
    odf$isgood <- FALSE

    # Ensure start position is at least 1
    odf$start <- pmax(odf$start, 1)

    #IF the overlap is somewhere NOT in the start of the read, we need the alignment to have mm < maxMismatch to be considered as good
    odf$isgood <- with(odf, ifelse(mismatch <= maxMismatch & start>1, TRUE, isgood))

    #now if the overlap is in the start of the read, we need mm<masMismatch AND the aligned portion has to be longer than length(marker)-1
    odf$isgood <- with(odf, ifelse(mismatch <= maxMismatch & start==1 & width>=nchar(marker)-1, TRUE, isgood))

    ##now for the good alignments, we trim from the start of the alignment
    #for the reads with no good alignment, we precautionarily trim half the length of the marker from the final of the read, to avoid having invomplete markers at the end
    #this is why we assign nchar(reads) - nchar(marker)/2 to all reads first, and then we change the value for the good alignments
    odf$cut <- nchar(reads) - floor(nchar(marker)/2)
    #we make sure that the cutting position is at least 1
    odf$cut <- with(odf, ifelse(isgood, pmax(odf$start-1, 1), cut))

    #we ensure the cut does not exceed read length AND that is not negative
    odf$cut <- pmin(odf$cut, nchar(reads))
    odf$cut <- pmax(odf$cut, 1)

    #now we convert the cut column into a vector to be able to trim the qual from there (((is this necessary?))
    cut_vector <- as.integer(odf$cut)

    # DEBUGGING checking invalid cut values, if there are some NA cut values, we stop the process
    if(any(is.na(cut_vector))) {
        stop("Invalid cut values detected (NA)")
    }

    #we trim the reads from the cut_vector point
    reads <- subseq(reads, 1, cut_vector)
    #we will convert the quality scores to a different format
    #first we convert it to a biological string to be able to apply the subseq function and get identical trimming from reads and qual
    qual_bstring <- BStringSet(as.character(qual@quality))
    qual_trimmed <- subseq(qual_bstring, 1, cut_vector)
    qual <- PhredQuality(as.character(qual_trimmed))

    #DEBUGGING
    if(length(reads) != length(qual)) {
        stop("Length mismatch inside function: reads=", length(reads), " qual=", length(qual))
    }
    if(any(width(reads) != width(qual))) {
        stop("Width mismatch inside function at positions: ", which(width(reads) != width(qual)))
    }

    return(list(reads = reads, qual = qual))

}

# Process both read files
result1 <- trim_overreading(reads1, qual1, largeLTRfrag_rc)
result2 <- trim_overreading(reads2, qual2, common_linker_rc)

trimmed_reads1 <- result1$reads
trimmed_qual1 <- result1$qual
trimmed_reads2 <- result2$reads
trimmed_qual2 <- result2$qual

#we only keep reads that are longer thatn 0 in both reads, to keep them paired
keep_indices <- (width(trimmed_reads1) > 0) & (width(trimmed_reads2) > 0)


# Create new ShortRead objects. we remove the ids of the reads that are empty
output_fq1 <- ShortReadQ(trimmed_reads1[keep_indices], trimmed_qual1[keep_indices], id(fq1)[keep_indices])
output_fq2 <- ShortReadQ(trimmed_reads2[keep_indices], trimmed_qual2[keep_indices], id(fq2)[keep_indices])

# Verify synchronization
if(length(output_fq1) != length(output_fq2)) {
  stop("Final synchronization failed: R1=", length(output_fq1), " R2=", length(output_fq2))
}

output_file1 <- paste0(meta, ".rc_removed_R1.fastq.gz")
output_file2 <- paste0(meta, ".rc_removed_R2.fastq.gz")

if (file.exists(output_file1)) file.remove(output_file1)
if (file.exists(output_file2)) file.remove(output_file2)

# Write output files using meta variable
writeFastq(output_fq1, output_file1, compress=TRUE)
writeFastq(output_fq2, output_file2, compress=TRUE)