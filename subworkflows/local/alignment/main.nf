include { ALIGNMENT_local } from '../../../modules/local/alignment/main'
include { GENOME_INDEXING_local } from '../../../modules/local/genome_indexing/main'

workflow ALIGNMENT_wfl {

    take:
    ch_dereplicated
    ch_refGenome             //take the genome name from the samplesheet, reference the samplesheet in the maindivided.nf

    main:
    // First, check which genomes need indexing
    ch_refGenome.branch { refGenome_name ->
        def index_path = "${params.runfolderDir}/../results/13_genome_index/${refGenome_name}_STAR_index"
        def index_exists = new File("${index_path}/SA").exists() && new File("${index_path}/SAindex").exists()
        
        log.info "Genome: ${refGenome_name}, Index exists: ${index_exists}"
        
        needs_indexing: !index_exists
            return refGenome_name
        has_index: index_exists
            return [refGenome_name, index_path]
    }.set { genome_branches }
    
    // Only run indexing for genomes that need it
    GENOME_INDEXING_local(genome_branches.needs_indexing)
    
    // Create genome index channel
    ch_existing_indexes = genome_branches.has_index
    ch_new_indexes = GENOME_INDEXING_local.out.index.map { index -> [refGenome_name, index] }
    ch_all_indexes = ch_existing_indexes.mix(ch_new_indexes)
    
    // Combine with dereplicated data
    ch_dereplicated
        .combine(ch_all_indexes.first())  // Use .first() to get only one genome index
        .map { sample, r1, r2, keys, refGenome_name, genome_index_dir ->
            return [sample, r1, r2, keys, genome_index_dir]
        }
        .set { ch_alignment_input }

    // Now pass the combined channel to alignment
    ALIGNMENT_local(ch_alignment_input)
    ALIGNMENT_local.out.aligned.set{ ch_aligned }
    // do we make different channels for different types of alignment? chimeras, multimapped, uniquely mapped?

    emit:
    ch_aligned

}