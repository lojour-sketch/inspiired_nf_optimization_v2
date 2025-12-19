include { ALIGNMENT_SINGLE_local } from '../../../modules/local/alignment/main_insp_until_alignment'
include { GENOME_INDEXING_local } from '../../../modules/local/genome_indexing/main'

workflow ALIGNMENT_wfl {

    take:
    ch_keys      //input from preprocessing workflow
    ch_refGenome             //take the genome name from the samplesheet, reference the samplesheet in the maindivided.nf

    main:

    //First we sort the keys channel by sample name
    ch_keys
        .toSortedList { a, b -> a[0] <=> b[0] }
        .flatMap { it }
        .set { ch_sorted_keys }

    //sorted_keys: tuple ( sample_name, r1, r2 ) in sample name order
    // Get unique genomes for indexing
    ch_unique_genomes = ch_refGenome
        .map { sample, refGenome_name, refGenomeFile, refKnowngeneFile -> 
            [refGenome_name, refGenomeFile] 
        }
        .unique()

    // Only run indexing for genomes that need it
    GENOME_INDEXING_local(ch_unique_genomes)

    // Combine with dereplicated data
    ch_sorted_keys
        .join ( ch_refGenome)
        .map { sample, r1, r2, keys, genome_name, refGenomeFile, refKnowngeneFile -> [genome_name, sample, r1, r2] }
        .set { ch_sample_with_genome }

    ch_sample_with_genome
        .combine( GENOME_INDEXING_local.out.index )  // we used combine{ ch_all_indexes, by: 0 } here in 25/10/23 and tried alignment with 144 samples
        .map { genome_name, sample, r1, r2, genome_name_duplicate, index_path -> [sample, r1, r2, index_path] }
        .set { ch_alignment_input}

    // Now pass the combined channel to alignment
    ALIGNMENT_SINGLE_local(ch_alignment_input)
    ALIGNMENT_SINGLE_local.out.aligned.set{ ch_aligned }
    // do we make different channels for different types of alignment? chimeras, multimapped, uniquely mapped?

    emit:
    ch_aligned

}