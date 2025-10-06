include { ALIGNMENT_local } from '../../../modules/local/alignment/main'
include { GENOME_INDEXING_local } from '../../../modules/local/genome_indexing/main'

workflow ALIGNMENT_wfl {

    take:
    ch_dereplicated
    genome_name             //take the genome name from the samplesheet, reference the samplesheet in the maindivided.nf

    main:
    // Check if index exists, only run indexing if needed, it takes quite long
    def index_path = "${params.runfolderDir}/../results/12_genome_index/${genome_name}_STAR_index"
    def index_exists = new File(index_path).exists()
    
    if (!index_exists) {
        GENOME_INDEXING_local(genome_name)
    }
    
    // Use existing index or newly created one
    // condition ? value_if_true : value_if_false
    def genome_index = index_exists ? index_path : GENOME_INDEXING_local.out."${genome_name}_STAR_index"

    ALIGNMENT_local(ch_dereplicated, genome_index)
    ALIGNMENT_local.out.aligned.set{ ch_aligned }
    // do we make different channels for different types of alignment? chimeras, multimapped, uniquely mapped?

    emit:
    ch_aligned

}