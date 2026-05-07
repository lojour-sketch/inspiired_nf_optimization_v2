process GENOME_INDEXING_local {

    publishDir "${params.outdir}/12_genome_index/${genome_name}", mode: 'copy', overwrite: true

    input:
    tuple val(genome_name), path(refGenomeFile)

    output:
    tuple val(genome_name), path("${genome_name}_STAR_index"), emit: index

     when:
    !file("${params.outdir}/${genome_name}/${genome_name}_STAR_index/SA").exists() || 
    !file("${params.outdir}/${genome_name}/${genome_name}_STAR_index/SAindex").exists()

    script:
    """
    mkdir -p ${genome_name}_STAR_index
    STAR --runMode genomeGenerate \
         --genomeDir ${genome_name}_STAR_index \
         --genomeFastaFiles ${params.runfolderDir}/../${refGenomeFile} \
         --runThreadN 8
    """
}