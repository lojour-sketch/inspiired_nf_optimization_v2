process GENOME_INDEXING_local {

    publishDir '/../../../results/12_genome_index', mode: 'symlink', overwrite: true

    input:
    val genome_name

    output:
    path "${genome_name}_STAR_index" dir true

    script:
    """
    module load STAR/2.7.11b-GCC-12.3.0
    mkdir -p ${genome_name}_STAR_index
    STAR --runMode genomeGenerate \
         --genomeDir ${genome_name}_STAR_index \
         --genomeFastaFiles ${genome_name}.fa \
         --runThreadN 8
    """
}