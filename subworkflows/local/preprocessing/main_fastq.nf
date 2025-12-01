include { DEMUXING_FASTQ_local } from '../../../modules/local/demuxing_fastq_local/main'
include { UMIEXTRACT_local } from '../../../modules/local/umi_extract_local/main'
include { FASTQCANDTRIM_wfl } from '../../../subworkflows/nf-core/fastqc_fastq_umitools_trimgalore/main'
include { MULTIQC_wfl } from '../../../subworkflows/nf-core/multiqc/main'
include { LTRchecking_seqkit_local } from '../../../modules/local/LTR_primer_checking_local/main_seqkit_fullfq'
include { RCremoval_inspiired_local } from '../../../modules/local/rcremoval_local/main_inspiired_like'
include { FINDVECTOR_local } from '../../../modules/local/find_vector_local/main'
include { SHORTREMOVE_local } from '../../../modules/local/short_seq_removal_local/main'

// MAIN

workflow PREPROCESSING_wfl {

    take:
    ch_linkers          //input type: [sample, unique_linker, common_linker]
    ch_primer_ltr       //input type: [sample, primer, ltrbit, largeLTRFrag, mingDNA]
    vectorfasta

    main:

    ///////////////////////// Passing undetermined fastq data into fastq per sample (demuxing) /////////////////////////
    ch_demux_fastq = Channel.of([params.FASTQfolderDir, file(params.demuxSampleSheet)])

    DEMUXING_FASTQ_local(ch_demux_fastq)
    DEMUXING_FASTQ_local.out.fastq.set { ch_demux_fastq }

    ///////////////////////// Extracting UMI from reads and adding it to headers /////////////////////////
        //crete necessary input channel for umi_extract_local
        //first create channel with tuple(sample_id, [R1, R2])
        ch_reads_by_sample = ch_demux_fastq
            .flatten()
            .map { file ->
                def name = file.getFileName().toString()
                def sample_name = name
                if (name.endsWith('.R1.fq.gz')) {
                    sample_name = name.substring(0, name.length() - '.R1.fq.gz'.length())
                } else if (name.endsWith('.R2.fq.gz')) {
                    sample_name = name.substring(0, name.length() - '.R2.fq.gz'.length())
                }
                tuple(sample_name, file)
            }
            .groupTuple(sort: true)

        // now we combine the reads channel with the linkers channel per sample
        ch_reads_by_sample
            .join(ch_linkers) //join by sample_id
            .map { sample_id, reads, linker1, linker2 ->
                    tuple(sample_id, linker1, linker2, reads)
            }
            .set { ch_umi_extract_input }
        //now running the umi_extract process
        UMIEXTRACT_local(ch_umi_extract_input)
            UMIEXTRACT_local.out.set { ch_umi_fastq}
        //changing umi output for FASTQCANDTRIM process. 
        def ch_umi_out = ch_umi_fastq
            .map { sample_id, fq1, fq2, log ->
                def meta = [ id: sample_id ]
                meta['single_end'] = (fq2 == null)
                tuple(meta, [fq1, fq2])
            }
        



    ///////////////////////// Filtering FASTQs by quality and removing Ns from sequences /////////////////////////
        //we will directly add the fastq_fastqc_umitools_trimgalore subworkflow
        FASTQCANDTRIM_wfl(
            ch_umi_out,           // input tuple: [ meta, reads ]
            false,                  //skip_fastqc
            false,                  // with_umi
            true,                   //skip_umi_extract
            false,                  // skip_trimming
            0,                   // umi_discard_read (0 , 1 or 2)
            10000,                   // min_trimmed_reads (default 10000)
        )

        def fastqc_html          = FASTQCANDTRIM_wfl.out.fastqc_html
        def fastqc_zip           = FASTQCANDTRIM_wfl.out.fastqc_zip
        def fastqc_html_trimmed  = FASTQCANDTRIM_wfl.out.fastqc_html_trimmed
        def fastqc_zip_trimmed   = FASTQCANDTRIM_wfl.out.fastqc_zip_trimmed
        def ch_fastq_filtered    = FASTQCANDTRIM_wfl.out.reads

        MULTIQC_wfl(
            fastqc_html,                 //rawfastqc html
            fastqc_zip,                 //rawfastqc zip
            fastqc_html_trimmed,                 //trimmedfastqc html
            fastqc_zip_trimmed,                //trimmedfastqc zip
            file('testresults/empty.yml'),                          //multiqc_config
            file('testresults/empty_extra.yml'),                          //extra_multiqc_config
            file('testresults/empty.png'),                          //multiqc_logo
            file('testresults/empty_replace.txt'),                          //replace_names
            file('testresults/empty_samples.txt'),                          //sample_names
        )




    ///////////////////////// Removing the LTR sequences from R2 and primers from R1 /////////////////////////
        //prepare input for LTRCHECKING_wfl
        ch_fastq_filtered
            .map { meta, reads -> tuple(meta.id, reads) }
            .join( ch_primer_ltr.map { id, primer, ltrbit, largeLTRFrag, project, mingDNA -> tuple(id, primer, ltrbit, largeLTRFrag, project, mingDNA) } )
            .map { sample_id, reads, primer, ltrbit, largeLTRFrag, project, mingDNA ->
                tuple(sample_id, reads, primer, ltrbit, largeLTRFrag, project, mingDNA)
            }
            .set { ch_input }
        
        LTRchecking_seqkit_local(ch_input)
        def ch_ltr_chunks = LTRchecking_seqkit_local.out.reads
    



    ///////////////////////// Removing the reverse complement reads, reverseLTRfrag form R1 and reversecommonlinker form R2 /////////////////////////
        //we need to first join the channels we need for rcremoval, if not nextflow will not cache properly
        ch_ltr_chunks
        .join( ch_linkers )
        .set { ch_joined_input }
        ch_ltr_chunks.view { "LTR chunks channel: ${it}" }

        RCremoval_inspiired_local(ch_joined_input)
        def ch_rc_removed = RCremoval_inspiired_local.out.reads
    



    ///////////////////////// Find vector sequences and remove them /////////////////////////
        ch_rc_removed
        .join( ch_primer_ltr )
        .set { ch_findvector_input }

        FINDVECTOR_local(ch_findvector_input, vectorfasta)
        def ch_vector_removed = FINDVECTOR_local.out.ch_vector_removed
    



    ///////////////////////// Remove short sequences /////////////////////////

        ch_vector_removed
        .join( ch_primer_ltr )
        .set { ch_shortremove_input }

        SHORTREMOVE_local(ch_shortremove_input)





    emit:
    ch_short_removed = SHORTREMOVE_local.out.reads


}