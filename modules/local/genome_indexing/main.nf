process GENOME_INDEXING_local {

    publishDir "${params.runfolderDir}/../results/13_genome_index", mode: 'symlink', overwrite: true

    input:
    tuple val(genome_name), path(refGenomeFile)

    output:
    tuple val(genome_name), path("${genome_name}_STAR_index"), emit: index

    script:
    """
    mkdir -p ${genome_name}_STAR_index
    STAR --runMode genomeGenerate \
         --genomeDir ${genome_name}_STAR_index \
         --genomeFastaFiles ${params.runfolderDir}/../${refGenomeFile} \
         --runThreadN 8
    """
}