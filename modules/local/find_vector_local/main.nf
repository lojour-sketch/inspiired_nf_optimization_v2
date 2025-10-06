process FINDVECTOR_local {

    conda "${params.runfolderDir}/../modules/local/find_vector/environment.yml"
    publishDir "${params.runfolderDir}/../results/10_findvector", mode: 'symlink', overwrite: true

    input:
    tuple val(meta), path(read1), path(read2)
    tuple val(meta), val(primer), val(ltrbit), val(largeLTRfrag), val(mingDNA)
    path(vector_fasta)

    output:
    tuple val(meta), path("${meta}.R1_vector_removed.fastq.gz"), path("${meta}.R2_vector_removed.fastq.gz"), emit: ch_vector_removed

    script:
    """
    Rscript -e '
    library(Biostrings)
    library(dplyr)
    
    vector <- "${vector_fasta}"
    read1 <- "${read1}"
    read2 <- "${read2}"
    primerltr <- "${primer}${ltrbit}"
    globalIdentity <- 0.75
    
    ## we will run the alignment with minimap. However, minimap calculates the identity percentage differently from BLAT, because minimap takes into account gaps to calculate the query length
    ## so we will calculate the identity the BLAT way (as inspiired) calculating matches / qSize from minimaps output
    ## then we will add this new identity to the output PAF file as the 13th column
    calculateidentity <- function(query_file, subject_file, min_identity) {
        # Run minimap2
        paf_file <- tempfile(pattern = "minimap", fileext = ".paf")
        system2("minimap2", 
                args = c("-x", "sr", "-N", "100", "-p", "0.1", "--eqx",
                        subject_file, query_file),
                stdout = paf_file, stderr = FALSE)
            
        if(file.exists(paf_file) && file.info(paf_file)\$size > 0) {
            paf <- read.table(paf_file, sep="\t", fill=TRUE, stringsAsFactors=FALSE)
            # we create a data frame that saves all the info
            result <- data.frame (
                qName = paf[,1],
                qSize = paf[,2], 
                qStart = paf[,3],
                matches = paf[,10],
                tStart = paf[,8],
                identity = paf[,10] / paf[,2],  # BLAT-style identity
                stringsAsFactors = FALSE
                )
            file.remove(paf_file)
            return(result)
        } else {
            return(data.frame())
        }
    }

    ## we align the reads
    hits.v.1 <- try(calculateidentity(read1, vector, globalIdentity)) 
    hits.v.2 <- try(calculateidentity(read2, vector, globalIdentity))
  
        
    ## we filter exactly as inspiired
    hits.v.1 <- dplyr::filter(hits.v.1, matches > globalIdentity*qSize & qStart <= 5)
    hits.v.2 <- dplyr::filter(hits.v.2, matches > globalIdentity*qSize)

    ## Extract common base name from headers (remove the " 1:N:0" / " 2:N:0" part)
    get_base_name <- function(header) {
        gsub(" [12]:N:0.*", "", header)
    }

    ## In the merge - use base names for matching
    hits.v.1\$baseName <- get_base_name(hits.v.1\$qName)
    hits.v.2\$baseName <- get_base_name(hits.v.2\$qName)

    ## now we merge the dataframes to get the reads that have vector in both pairs. we will only remove these    
    hits.v <- try(merge(hits.v.1[, c("baseName", "tStart")],
                        hits.v.2[, c("baseName", "tStart")],
                        by="baseName")
                    ,silent = TRUE)
    ## as multihits can happen, we take only unique names
    vqName_base <- unique(hits.v\$baseName) 
        
    #now we remove reads from both files. first we read them and get their basenames
    reads1 <- readDNAStringSet(read1)
    reads2 <- readDNAStringSet(read2)

    reads.1_base <- get_base_name(names(reads1))
    reads.2_base <- get_base_name(names(reads2))

    ## Remove vector reads from both files
    ## vqName_base contains the base names of the headers that belong to the pairs that both have vectors 
    reads.1_clean <- reads1[!reads.1_base %in% vqName_base]
    reads.2_clean <- reads2[!reads.2_base %in% vqName_base]

    writeXStringSet(reads.1_clean, '${meta}.R1_vector_removed.fastq.gz', format='fastq', compress=TRUE)
    writeXStringSet(reads.2_clean, '${meta}.R2_vector_removed.fastq.gz', format='fastq', compress=TRUE)

    '

    """

}