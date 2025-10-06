process DEREPLICATE_local {

    publishDir "${params.runfolderDir}/../results/12_dereplicated_reads", mode: 'copy', overwrite: true

    input:
    tuple val(sample), path(read1), path(read2)

    output:
    tuple val(sample), path("${sample}.R1_unique_by_sequence.fastq.gz"), path("${sample}.R2_unique_by_sequence.fastq.gz"), path("${sample}_keys.RData"), emit: ch_dereplicated

    script:
    """
    Rscript -e 'library(Biostrings)
      r1 <- readDNAStringSet("${read1}", format="fastq")
      r2 <- readDNAStringSet("${read2}", format="fastq")
      length(r1)
      length(r2)
      r1u <- unique(r1)
      r2u <- unique(r2)
      length(r1u)
      length(r2u)
      writeXStringSet(r1u, "${sample}.R1_unique_by_sequence.fastq", format="fastq")
      writeXStringSet(r2u, "${sample}.R2_unique_by_sequence.fastq", format="fastq")
      # we create a dataframe that will contain the read name, the reads unique r1 match and unique r2 match (the number/position of the unique reads in the set)
      keys <- data.frame(name=names(r1),
                         R1=match(as.character(r1), as.character(r1u)), # match finds the index of the unique sequence in r1u that matches each original read in r1. 
                         R2=match(as.character(r2), as.character(r2u)))
      #taking the corresponding unique reads for each read, or more precisely the index of the unique reads, we form a paired key like 3_5 
      #(read1 is the 3rd unique read, read2 is the 5th uique read)
      keys\$readPairKey <- paste0(keys\$R1, "_", keys\$R2)
      save(keys, file="${sample}_keys.RData")'

      #compressing the output files
      gzip ${sample}.R1_unique_by_sequence.fastq
      gzip ${sample}.R2_unique_by_sequence.fastq
    """

}