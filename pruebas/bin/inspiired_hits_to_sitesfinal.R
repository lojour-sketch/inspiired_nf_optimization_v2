#!/usr/bin/env Rscript


#load libraries
library(GenomicRanges)
library(IRanges)

#load functions
dereplicateSites <- function(uniqueReads){
    if(length(uniqueReads) == 0) return(print("No sites to dereplicate"))
    #do the dereplication, but loose the coordinates
    sites.reduced <- reduce(flank(uniqueReads, -5, both=TRUE), with.revmap=T)
    sites.reduced$counts <- sapply(sites.reduced$revmap, length)
    
    #order the unique sites as described by revmap
    dereplicatedSites <- uniqueReads[unlist(sites.reduced$revmap)]
    
    #if no sites are present, skip this step - keep doing the rest to provide a
    #similar output to a successful dereplication
    if(length(uniqueReads)>0){
      #split the unique sites as described by revmap (sites.reduced$counts came from revmap above)
      dereplicatedSites <- split(dereplicatedSites, Rle(values=seq(length(sites.reduced)), lengths= sites.reduced$counts))
    }
    
    #do the standardization - this will pick a single starting position and
    #choose the longest fragment as ending position
    dereplicatedSites <- unlist(reduce(dereplicatedSites, min.gapwidth=5))
    mcols(dereplicatedSites) <- mcols(sites.reduced)
    
    dereplicatedSites
  }


#load arguments
args <- commandArgs(trailingOnly = TRUE)
hits.r1 <- args[1] #tsv
hits.r2 <- args[2] #tsv
sitesfinal_file <- args[3] #RData



sample <- args[4] #character
keys_file <- args[5] #RData

sample2 <- gsub("_", "-", sample)

#load keys
object <- load(keys_file)
keys <- get(object)

#hits are already filtered byquality: maxAlignstart, minPctIdent, maxFragLength
hits.R1_tsv <- read.table(hits.r1, header=TRUE, sep="\t", stringsAsFactors=FALSE, row.names=1, colClasses=c("character", "character", "integer", "integer", "integer", "character", "integer", "integer", "integer", "integer", "character", "integer", "character"))
hits.R2_tsv <- read.table(hits.r2, header=TRUE, sep="\t", stringsAsFactors=FALSE, row.names=1, colClasses=c("character", "character", "integer", "integer", "integer", "character", "integer", "integer", "integer", "integer", "character","integer", "character"))

message("Creating GRanges from hits.R1_tsv and hits.R2_tsv")

# Rename conflicting columns in the data frames BEFORE creating GRanges, to avoid duplicate column name errors
reserved <- c("seqnames","ranges","strand","seqlevels","seqlengths","isCircular","start","end","width","element")

rename_conflicting_cols <- function(df, reserved) {
  conflicting <- intersect(names(df), reserved)
  for (col in conflicting) {
    names(df)[names(df) == col] <- paste0("orig_", col)
  }
  return(df)
}

hits.R1_tsv <- rename_conflicting_cols(hits.R1_tsv, reserved)
hits.R2_tsv <- rename_conflicting_cols(hits.R2_tsv, reserved)

hits.R1 <- makeGRangesFromDataFrame(hits.R1_tsv, 
                                   keep.extra.columns = TRUE,
                                   seqnames.field = "orig_seqnames",
                                   start.field = "orig_start", 
                                   end.field = "orig_end",
                                   strand.field = "orig_strand",
                                   starts.in.df.are.0based = FALSE)

message("Created GRanges for hits.R1")
message(paste0("Number of hits.R1: ", length(hits.R1)))

hits.R2 <- makeGRangesFromDataFrame(hits.R2_tsv, 
                                   keep.extra.columns = TRUE,
                                   seqnames.field = "orig_seqnames",
                                   start.field = "orig_start", 
                                   end.field = "orig_end",
                                   strand.field = "orig_strand",
                                   starts.in.df.are.0based = FALSE)

message("Created GRanges for hits.R2")
message(paste0("Number of hits.R2: ", length(hits.R2)))

#we will identifyall combinations of unique r1 and r2 sequences present in the data

unique_key_pairs <- unique(keys[,c("R1", "R2", "readPairKey")])

message("Identified unique key pairs")
message(paste0("Number of unique key pairs: ", nrow(unique_key_pairs)))
#now we will reduce the alignments
#this way, we will identify the distinct genomic locations `resent in the data for the 
# R1 sequences (breakpoint positions)
# R2 sequences (integration site positions)
#we get the unique locis
red.hits.R1 <- reduce(flank(hits.R1, -1, start=TRUE), min.gapwidth=0L, with.revmap=TRUE)
red.hits.R2 <- reduce(flank(hits.R2, -1, start=TRUE), min.gapwidth=0L, with.revmap=TRUE)

message("Reduced hits.R1 and hits.R2 to unique loci")
message(paste0("Number of reduced hits.R1: ", length(red.hits.R1)))
message(paste0("Number of reduced hits.R2: ", length(red.hits.R2)))
#now we will find all possible combinations of R1 and R2 which meet criteria for pairing
#    each pairing must come form one R1 and one R2
#    paired loci should be present in opposite strands
#    correct downstream orientation

    # first we use findOverlaps to find the overlaps and check the first criteria
    pairs <- findOverlaps(
        red.hits.R1,
        red.hits.R2,
        maxgap=2500,
        ignore.strand=TRUE
    )

    message("Found overlaps between reduced hits.R1 and hits.R2")
    # we will use the overlaps to find the R1 and R2 loci
    R1.loci <- red.hits.R1[queryHits(pairs)]
    R2.loci <- red.hits.R2[subjectHits(pairs)]

    # now we check if the R1 and R2 are in opposite strands and in correct orientation
    R1.loci.starts <- start(R1.loci)
    R2.loci.starts <- start(R2.loci)
    R1.loci.strand <- strand(R1.loci)
    R2.loci.strand <- strand(R2.loci)

    keep.loci <- ifelse(
        R2.loci.strand == "+",
        as.vector(R1.loci.starts > R2.loci.starts &
                    R1.loci.strand != R2.loci.strand),
        as.vector(R1.loci.starts < R2.loci.starts &
                    R1.loci.strand != R2.loci.strand)
    )

    keep.loci <- as.vector(
        keep.loci & R2.loci.strand!= "*" & R1.loci.strand != "*"
    )

    R1.loci <- R1.loci[keep.loci]
    R2.loci <- R2.loci[keep.loci]

message("Identified valid loci that meet pairing criteria")
message(paste0("Number of valid R1 locis: ", length(R1.loci)))
message(paste0("Number of valid R2 locis: ", length(R2.loci)))

#now that we have the correct locis, we will create a genomic loci key which links genomic loci to R1 and R2 sequences
    loci.key <- data.frame(
        "R1.loci" = queryHits(pairs)[keep.loci],
        "R2.loci" = subjectHits(pairs)[keep.loci]
    )
    #we add the pairkeys that the locis of r1 and r2 create
    loci.key$lociPairKey <- paste0(loci.key$R1.loci, ":", loci.key$R2.loci)

    #we now add the query names of r1 and r2 sequences
    #this means, the query names of the reads with those loci
    loci.key$R1.qNames <- IntegerList(lapply(R1.loci$revmap, function(x){
        as.integer(mcols(hits.R1)$qName[x])
    }))

    loci.key$R2.qNames <- IntegerList(lapply(R2.loci$revmap, function(x){
        as.integer(mcols(hits.R2)$qName[x])
    }))

    #we now add the readpairs that contain each r1 or r2 query name
    #this is, for each r1 or r2 read we add the readpairs that contain the querynames we already added (the ones with the corresponding loci)

    message("debugging")
    message("A part of loci.key$R1.qNames (first 10 elements):")
    print(head(loci.key$R1.qNames, 10))

    message("A part of unique_key_pairs$R1 (first 10):")
    print(head(unique_key_pairs$R1, 10))

    message("A part of loci.key$R2.qNames (first 10):")
    print(head(loci.key$R2.qNames, 10))

    message("A part of unique_key_pairs$R2 (first 10):")
    print(head(unique_key_pairs$R2, 10))  
    loci.key$R1.readPairs <- IntegerList(lapply(
        loci.key$R1.qNames, function(x){
            which(unique_key_pairs$R1 %in% x)
    }))

    loci.key$R2.readPairs <- IntegerList(lapply(
        loci.key$R2.qNames, function(x){
            which(unique_key_pairs$R2 %in% x)
    }))

    message("Created loci key linking R1 and R2 loci to read pairs")
    message(paste0("Number of loci pairs: ", nrow(loci.key)))

    message("debugging")
    message(paste0("R1 qNames in loci.key: ", length(loci.key$R1.qNames)))
    message(paste0("Names in unique_key_pairs: ", length(unique_key_pairs$R1)))
    message(paste0("R2 qNames in loci.key: ", length(loci.key$R2.qNames)))
    message(paste0("Names in unique_key_pairs: ", length(unique_key_pairs$R2)))
#we will now create a GRanges object from R1 and R2
#R2 contains the insertion sites
#R1 contains the various breakpoints

    #we create granges object with united r1 and r2 loci, with the loci pair key
    paired.loci <- GRanges(
        seqnames = seqnames(R2.loci),
        ranges = IRanges(
            start = ifelse(strand(R2.loci) == "+", start(R2.loci), start(R1.loci)),
            end = ifelse(strand(R2.loci) == "+", end(R1.loci), end(R2.loci))),
        strand = strand(R2.loci)
    )

    message(paste0("Number of paired loci before adding ReadPairKeys: ", length(paired.loci)))

    mcols(paired.loci)$lociPairKey <- loci.key$lociPairKey
    #we add the readpair keys that the read pairs that are present in unique_key_pairs contain
    #ONLY the pairs present in unique_key_pairs
    mcols(paired.loci)$readPairKeys <- CharacterList(lapply(
        1:length(paired.loci), 
        function(i){
            unique_key_pairs[intersect(
            loci.key$R1.readPairs[[i]], 
            loci.key$R2.readPairs[[i]]),
            "readPairKey"]
    }))

    # We emove R1:R2 pairings that do not have a readpairkey and are not present in the sequence data
    message("debugging")
    message("This is a part of loci.key$R1.readPairs")
    print(loci.key$R1.readPairs[1:5])
    message("This is a part of loci.key$R2.readPairs")
    print(loci.key$R2.readPairs[1:5])
    message("Filtering paired.loci to only those with intersected readPairKeys between R1 and R2")
    paired.loci <- paired.loci[sapply(paired.loci$readPairKeys, length) > 0]

    message("created paired.loci, linking loci pairs to read pairs")
    message(paste0("Number of paired loci after adding metadata: ", length(paired.loci)))

    #now we expand by readPairKeys. We will get a dataframe with all the readpairkeys (expanded, not in a list) and their corresponding locis (which can be repeated)
    ######changed from inspiired because length parameter of Rle was giving errors
    read.loci.mat <- data.frame(
        "lociPairKey" = rep(paired.loci$lociPairKey, sapply(paired.loci$readPairKeys, length)),
        "readPairKey" = unlist(paired.loci$readPairKeys)
    )

#we will only use the alignments that are termed unique
    # for that first we count all the pairs and how many times they appear
    readPairCounts <- table(read.loci.mat$readPairKey)
    #if the pairs appear more than once, we will take them as multihits or chimeras, so we will remove them
    uniq.readPairs <- names(readPairCounts[readPairCounts == 1])

    message("Identified unique read pairs")
    message(paste0("Number of unique read pairs: ", length(uniq.readPairs)))

############################# INSPIIRED NOW DETECTS CHIMERAS AND MULTIHITS, WE WILL ONLY DETECT UNIQUE SITES ##############################
    # from the matrix of all the readpairs, we will onlt keep the ones that appear in the uniq.readPairs list
    uniq.read.loci.mat <- read.loci.mat[read.loci.mat$readPairKey %in% uniq.readPairs,]
    message("Identified unique read loci")
    message(paste0("Number of unique read loci: ", length(uniq.read.loci.mat$readPairKey)))
    # take from paired.loci (united r1r2 with paired loci of target) the pairs with locis that only appear once (present in uniq.read.loci.mat)
    uniq.templates <- paired.loci[match(uniq.read.loci.mat$lociPairKey, paired.loci$lociPairKey)]
    uniq.templates$readPairKeys <- NULL
    #take the readpairs from those unique loci
    uniq.templates$readPairKey <- uniq.read.loci.mat$readPairKey
    message("Identified unique templates")
    message(paste0("Number of unique templates: ", length(uniq.templates)))

    # only take the unique keys (number mixes)
    uniq.keys <- keys[keys$readPairKey %in% uniq.readPairs,]
    # in uniq.reads we have the pairs that contain the keys that appear both in uniq.templates and uniq.keys
    #then we will add the samplename and ID of the reads
    uniq.reads <- uniq.templates[
        match(uniq.keys$readPairKey, uniq.templates$readPairKey)
    ]
    message("Identified unique keys")
    message(paste0("Number of unique keys: ", length(uniq.keys)))
    names(uniq.reads) <- as.character(uniq.keys$names)
    uniq.reads$sampleName <- sapply(
        strsplit(as.character(uniq.keys$names), "%"), "[[", 1
    )
    uniq.reads$ID <- sapply(strsplit(as.character(uniq.keys$names), "%"), "[[", 2)

    allSites <- uniq.reads
    message("Identified all sites")
    message(paste0("Number of all sites: ", length(allSites)))

    save(allSites, file=paste0(sample, "_allSites.RData"))

    sites.final <- dereplicateSites(allSites)
    message("Identified sites.final")
    message(paste0("Number of sites.final: ", length(sites.final)))
    if(length(sites.final)>0){
        sites.final$sampleName <- allSites[1]$sampleName
        sites.final$posid <- paste0(as.character(seqnames(sites.final)),
                                    as.character(strand(sites.final)),
                                    start(flank(sites.final, width=-1, start=TRUE)))
    }

    message(paste0("Number of sites.final after adding posid: ", length(sites.final)))
    message("Saving sites.final to RData")

    save(sites.final, file=paste0(sample, "_sitesfinal.RData"))