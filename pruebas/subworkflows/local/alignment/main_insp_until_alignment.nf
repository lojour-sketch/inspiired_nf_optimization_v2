include { ALIGNMENT_SINGLE_local } from '../../../modules/local/alignment/main_insp_until_alignment'
include { GENOME_INDEXING_local } from '../../../modules/local/genome_indexing/main'

workflow ALIGNMENT_wfl {

    take:
    ch_refGenome             //take the genome name from the samplesheet, reference the samplesheet in the maindivided.nf

    main:
    // Get unique genomes for indexing
    ch_unique_genomes = ch_refGenome
        .map { sample, refGenome_name, refGenomeFile, refKnowngeneFile -> 
            [refGenome_name, refGenomeFile] 
        }
        .unique()

    // Only run indexing for genomes that need it
    GENOME_INDEXING_local(ch_unique_genomes)

    Channel
        .fromPath("${params.runfolderDir}/../demoDataSet/*/", type: 'dir')
        .filter { it.isDirectory() }
        .filter { folder -> 
            folder.getName() ==~ /GTSP\d+-\d+/  // ← Only match GTSP####-# pattern
        }
        .map { folder ->
            def sample = folder.getName()
            def sampleFilePrefix = sample.replaceAll('-', '_')
            def r1_file = file("${folder}/${sampleFilePrefix}.R1.fa")
            def r2_file = file("${folder}/${sampleFilePrefix}.R2.fa")
            if (!r1_file.exists() || !r2_file.exists()) {
                log.warn "Missing merged files for sample ${sample}: R1=${r1_file.exists()}, R2=${r2_file.exists()}"
                return null
            }
            [sampleFilePrefix, r1_file, r2_file]
        }
        .filter { it != null }
        .set { ch_samples_with_reads }
    
    // Combine with dereplicated data
    ch_samples_with_reads
        .join ( ch_refGenome, by: 0 )
        .map { sample, r1, r2, genome_name, refGenomeFile, refKnowngeneFile -> [genome_name, sample, r1, r2] }
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