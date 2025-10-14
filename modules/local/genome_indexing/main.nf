process GENOME_INDEXING_local {

    publishDir '/../../../results/13_genome_index', mode: 'symlink', overwrite: true

    input:
    val genome_name

    output:
    path "${genome_name}_STAR_index", emit: index

    script:
    """
    mkdir -p ${genome_name}_STAR_index
    STAR --runMode genomeGenerate \
         --genomeDir ${genome_name}_STAR_index \
         --genomeFastaFiles ${params.runfolderDir}/../${genome_name}.fa \
         --runThreadN 8
    """
}