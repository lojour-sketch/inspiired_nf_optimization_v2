process RCremoval_inspiired {

    cpus 12
    cache 'deep'
    memory '50GB'

    publishDir '/home/lrenteria/inspiired_nf/results/9_reverse_complement_removal_inspiired', mode: 'symlink', overwrite: true

    input:
    tuple val(meta), path(read1), path(read2), val(primer), val(ltrbit), val(largeLTRfrag), val(mingDNA)
    tuple val(sample), val(unique_linker), val(common_linker)

    output:
    tuple val(meta), path("${meta}.rc_removed_R1.fastq.gz"), path("${meta}.rc_removed_R2.fastq.gz"), emit: reads
    tuple val(meta), val(primer), val(ltrbit), val(largeLTRfrag), val(mingDNA)

    script:
    """
    Rscript -e 'library(Biostrings)
    library(ShortRead)
    library(pwalign) #we need it for the pairwiseAlignment and nucleotideSubstitutionMatrix functions

    fq1 <- readFastq("${read1}")
    fq2 <- readFastq("${read2}")

    #R needs to convert the quality and the sequence to an specific format
    reads1 <- sread(fq1)
    reads2 <- sread(fq2)

    qual1 <- quality(fq1)
    qual2 <- quality(fq2)

    #we obtain the reverse complements of the markers
    largeLTRfrag <- "${largeLTRfrag}"
    common_linker <- "${common_linker}"
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
                                 type="overlap",)
        
        #we convert the pairwisealignment output into a dataframe
        odf <- PairwiseAlignmentsSingleSubject2DF(tmp)

        #we assume no reads have good alignments
        odf\$isgood <- FALSE

        #IF the overlap is somewhere NOT in the start of the read, we need the alignment to have mm < maxMismatch to be considered as good
        odf\$isgood <- with(odf, ifelse(mismatch <= maxMismatch & start>1, TRUE, isgood))

        #now if the overlap is in the start of the read, we need mm<masMismatch AND the aligned portion has to be longer than length(marker)-1
        odf\$isgood <- with(odf, ifelse(mismatch <= maxMismatch & start==1 & width>=nchar(marker)-1, TRUE, isgood))

        ##now for the good alignments, we trim from the start of the alignment
        #for the reads with no good alignment, we precautionarily trim half the length of the marker from the final of the read, to avoid having invomplete markers at the end
        #this is why we assign nchar(reads) - nchar(marker)/2 to all reads first, and then we change the value for the good alignments
        odf\$cut <- nchar(reads) - floor(nchar(marker)/2)
        odf\$cut <- with(odf, ifelse(isgood, odf\$start-1, cut))

        #now we convert the cut column into a vector to be able to trim the qual from there (((is this necessary?))
        cut_vector <- as.integer(odf\$cut)

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

    trimmed_reads1 <- result1\$reads
    trimmed_qual1 <- result1\$qual
    trimmed_reads2 <- result2\$reads
    trimmed_qual2 <- result2\$qual



    # Create new ShortRead objects
    output_fq1 <- ShortReadQ(trimmed_reads1, trimmed_qual1, id(fq1))
    output_fq2 <- ShortReadQ(trimmed_reads2, trimmed_qual2, id(fq2))

    # we will remove any possible output files form other runs
    if (file.exists("${meta}.rc_removed_R1.fastq.gz")) file.remove("${meta}.rc_removed_R1.fastq.gz")
    if (file.exists("${meta}.rc_removed_R2.fastq.gz")) file.remove("${meta}.rc_removed_R2.fastq.gz")

    # Write output files
    writeFastq(output_fq1, "${meta}.rc_removed_R1.fastq.gz", compress=TRUE)
    writeFastq(output_fq2, "${meta}.rc_removed_R2.fastq.gz", compress=TRUE)'

    """
}