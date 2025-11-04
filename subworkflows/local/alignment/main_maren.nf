include { ALIGNMENT_local } from '../../../modules/local/alignment/main'
include { GENOME_INDEXING_local } from '../../../modules/local/genome_indexing/main_maren'

workflow ALIGNMENT_wfl {

    take:
    ch_short_removed
    ch_refGenome             //take the genome name from the samplesheet, reference the samplesheet in the maindivided.nf

    main:
    // First, check which genomes need indexing
    ch_refGenome.branch { refGenome_name, refGenomeFile, refKnowngeneFile->
        def index_path = "${params.runfolderDir}/../results/13_genome_index/${refGenome_name}_STAR_index"
        def index_exists = new File("${index_path}/SA").exists() && new File("${index_path}/SAindex").exists()
        
        
        needs_indexing: !index_exists
            return [refGenome_name, refGenomeFile]
        has_index: index_exists
            return [refGenome_name, index_path]
    }.set { genome_branches }

    // Get unique genomes that need indexing, we only need one index per genome
    ch_unique_genomes_to_index = genome_branches.needs_indexing.unique()


    // Only run indexing for genomes that need it
    GENOME_INDEXING_local(ch_unique_genomes_to_index)
    
    // Create genome index channel
    ch_existing_indexes = genome_branches.has_index
    ch_new_indexes = GENOME_INDEXING_local.out.index
    ch_all_indexes = ch_existing_indexes.mix(ch_new_indexes)

    // Combine with dereplicated data
    ch_short_removed
        .merge ( ch_refGenome )
        .map { sample, r1, r2, genome_name, refGenomeFile, refKnowngeneFile -> [genome_name, sample, r1, r2] }
        .set { ch_sample_with_genome }

    ch_sample_with_genome
        .join( ch_all_indexes )  // we used combine{ ch_all_indexes, by: 0 } here in 25/10/23 and tried alignment with 144 samples
        .map { genome_name, sample, r1, r2, index_path -> [sample, r1, r2, index_path] }
        .set { ch_alignment_input}


    // Now pass the combined channel to alignment
    ALIGNMENT_local(ch_alignment_input)
    ALIGNMENT_local.out.aligned.set{ ch_aligned }
    // do we make different channels for different types of alignment? chimeras, multimapped, uniquely mapped?

    emit:
    ch_aligned

}