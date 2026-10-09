args <- commandArgs(trailingOnly=TRUE)
dir.create(args[1],recursive=TRUE,showWarnings=FALSE);setwd(args[1])
root <- normalizePath(args[2])
.libPaths(c(file.path(root,"vendor/bushman_r_library"),.libPaths()))
suppressPackageStartupMessages({library(Biostrings);library(ShortRead);library(rtracklayer)})
set.seed(197)
genome <- paste(sample(c("A","C","G","T"),24000,replace=TRUE),collapse="")
substr(genome,14001,15000) <- substr(genome,10001,11000)
dna <- DNAStringSet(c(chrSynthetic=genome))
export(dna,"reference.2bit")
vector <- DNAStringSet(c(vector=paste(sample(c("A","C","G","T"),500,replace=TRUE),collapse="")))
writeXStringSet(vector,"vector.fa")
primer <- "GAAAATC";ltr <- "TCTAGCA";linker <- "CTCCGCTTAAGGGACT"
rc <- function(x) as.character(reverseComplement(DNAString(x)))
sheet <- list();manifest <- list()
for(sid in c("A","B","Empty")) {
    r1 <- r2 <- character();names <- character()
    if(sid!="Empty") for(junction in c(2001L,5001L,8001L,10001L,18000L)) {
        for(i in seq_len(24L)) {
            len <- sample(220:450,1)
            opposite <- if(junction==18000L) junction-len+1L else junction+len-1L
            g1 <- if(junction==18000L) substr(genome,opposite,opposite+89L) else rc(substr(genome,opposite-89L,opposite))
            g2 <- if(junction==18000L) rc(substr(genome,junction-89L,junction)) else substr(genome,junction,junction+89L)
            r1 <- c(r1,paste0("ACGTACGTACGT",linker,g1))
            r2 <- c(r2,paste0(primer,ltr,g2))
            names <- c(names,paste0(sid,"_read_",length(names)+1L))
        }
    }
    if(sid!="Empty") for(i in seq_len(3L)) {
        r1 <- c(r1,paste0("ACGTACGTACGT",linker,rc(substr(genome,20911L,21000L))))
        r2 <- c(r2,paste0(primer,ltr,substr(genome,2001L,2090L)))
        names <- c(names,paste0(sid,"_chimera_",i))
    }
    for(mate in c("R1","R2")) {
        reads <- if(mate=="R1") r1 else r2
        fq <- ShortReadQ(DNAStringSet(reads),PhredQuality(BStringSet(vapply(nchar(reads),
            function(n) paste(rep("I",n),collapse=""),character(1)))),BStringSet(names))
        target <- paste0(sid,".",mate,".fq.gz")
        if(file.exists(target)) unlink(target)
        writeFastq(fq,target,compress=TRUE)
    }
    sheet[[sid]] <- data.frame(Sample_ID=sid,index="AAAAAAAAAAAAAAAAAAAA",
        index2=paste(rep(if(sid=="A") "A" else if(sid=="B") "C" else "G",12),collapse=""),
        common_linker=linker,primer=primer,ltrbit=ltr,
        largeLTRFrag="TGCTAGAGATTTTCCACACTGACTAAAAGGGTCTG",Sample_Project="truth",
        mingDNA=30,minPctIdent=97,maxAlignStart=5,maxFragLength=2500,refGenome="hg38",
        vectorSeq=normalizePath("vector.fa"),Specimen_ID=if(sid=="Empty") "Empty" else "S1",Replicate_ID=sid)
    manifest[[sid]] <- data.frame(Sample_ID=sid,r1=normalizePath(paste0(sid,".R1.fq.gz")),
        r2=normalizePath(paste0(sid,".R2.fq.gz")))
}
write.csv(do.call(rbind,sheet),"samples.csv",row.names=FALSE)
write.csv(do.call(rbind,manifest),"fastqs.csv",row.names=FALSE)
jsonlite::write_json(list(hg38=list(twobit=normalizePath("reference.2bit"),
    sequence_policy="synthetic truth: one repeated segment in the same reference")),"references.json",pretty=TRUE,auto_unbox=TRUE)
