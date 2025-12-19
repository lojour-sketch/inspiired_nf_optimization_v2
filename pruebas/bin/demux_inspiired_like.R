#! /usr/bin/env Rscript

library(ShortRead)

#load arguments
args <- commandArgs(trailingOnly=TRUE)
FASTQfolder <- args[1]
SampleSheet <- args[2]
indices <- strsplit(args[3], ",")[[1]]
samples <- strsplit(args[4], ",")[[1]]

print(samples)
print(indices)

#load functions
demultiplex_reads <- function(reads, suffix, I1Names, samples, indices) {
    I1Seqs <- as.vector(sread(I1))
    message("See samples: ", samples)
    message("See indices: ", indices)
    for (j in seq_along(samples)){
        sample_name <- samples[j]
        barcode <- indices[j]

        message(paste0('Demultiplexing ', suffix, ' read: ', j, '/', length(samples)))

        reads_for_sample <- reads[I1Seqs == barcode]
        fq_file <- paste0(sample_name, "_", suffix, ".fastq.gz")
        if (length(reads_for_sample) > 0) {
            writeFastq(reads_for_sample, fq_file, mode="w")
        } else {
            message("Warning: no reads found for sample", sample_name)
        }
    }  
}

    I1_files <- list.files(FASTQfolder, pattern="correctedI1-.", full.names=TRUE)
    I1 <- readFasta(I1_files)

  
  I1 <- I1[as.vector(sread(I1)) %in% indices]
  
   #only necessary if using native data - can parse out description w/ python
   I1Names <-  sapply(strsplit(as.character(ShortRead::id(I1)), " "), "[[", 1)#for some reason we can't dynamically set name/id on ShortRead!


message("See samples before demultiplexing: ", samples)
message("See indices before demultiplexing: ", indices)

  message('Starting to demultiplex R1')
  R1 <- readFastq(paste0(FASTQfolder,"/Undetermined_S0_L001_R1_001.fastq.gz"))
  demultiplex_reads(R1, "R1", I1Names, samples, indices)
  message('completed demultiplexing R1')  

  message('Starting to demultiplex R2')
  R2 <- readFastq(paste0(FASTQfolder,"/Undetermined_S0_L001_R2_001.fastq.gz"))
  demultiplex_reads(R2, "R2", I1Names, samples, indices)
  message('completed demultiplexing R2')
