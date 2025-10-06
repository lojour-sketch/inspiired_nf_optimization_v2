#!/usr/bin/env nextflow

// Parameters
params.samplesheet = '' // this sheet contains barcodes, processing parameters, LTR fragments, etc.
params.runfolderDir = ''
params.outputDir = ''
params.linkerdata = ''

// Parameter validation
//if (!params.samplesheet) error "Missing --samplesheet parameter, provide path to samplesheet"
if (!params.runfolderDir) error "Missing --runfolderDir parameter. provide path to BCL data folder or tarball"
if (!params.outputDir) error "Missing --outputDir parameter, provide path to output directory"
if (!params.linkerdata) error "Missing --linkerdata parameter, provide .tsv file with linker sequence information"

// Log info

log.info "Using data from directory: ${params.runfolderDir}"

// Include modules

include { BCL2FASTQ_wfl } from './subworkflows/local/bcl2fastq/main.nf'
include { EXTRACTUMI_wfl } from './subworkflows/local/umi_extract/main.nf'
include { FASTQCANDTRIM_wfl } from './subworkflows/nf-core/fastq_fastqc_umitools_trimgalore/main.nf'
include { MULTIQC_wfl } from './subworkflows/local/multiqc/main.nf'
include { LTRCHECKING_wfl } from './subworkflows/local/LTR_primer_checking/main' 
include { RCREMOVAL_wfl } from './subworkflows/local/rcremoval/main'
include { SHORTREMOVE_wfl } from './subworkflows/local/short_seq_removal/main'
include { DEREPLICATE_wfl } from './subworkflows/local/dereplicate/main'
include { FINDVECTOR_wfl } from './subworkflows/local/findvector/main'

// Create necessary input tuple for bcl2fastq

Channel
    .of( tuple([id: 'run1'], file(params.samplesheet), file(params.runfolderDir, type: 'dir')) )
    .set { ch_bcl_input }



// Workflow
workflow {
    // ******************************** NECESSARY CHANNELS ********************************
        // we need a channel with the sample unique and common linkers to extract the umi
            Channel
            .fromPath(params.samplesheet)
            .map { file ->
                def lines = file.text.readLines()
                def dataStart = lines.findIndexOf { it.startsWith('[Data]') } + 1
                def csvLines = lines[dataStart..-1].join("\n")
                return csvLines
            }
            .splitCsv(header: true)
            .map { row ->
                tuple(row.Sample_ID, row.sample_unique_linker, row.common_linker)
            }
            .set { ch_linkers }
        // we need a channel with some samplesheet information as input for some of the processes
        Channel
            .fromPath(params.samplesheet)
            .map { file ->
                def lines = file.text.readLines()
                def dataStart = lines.findIndexOf { it.startsWith('[Data]') } + 1
                def csvLines = lines[dataStart..-1].join("\n")
                return csvLines
            }
            .splitCsv(header: true)
            .map { row ->
                tuple(row.Sample_ID, row.primer, row.ltrbit, row.largeLTRFrag, row.mingDNA)
            }
            .set { ch_primer_ltr }

    // ************************************* WORKFLOW *************************************
    BCL2FASTQ_wfl(ch_bcl_input)
    def ch_demux_fastq = BCL2FASTQ_wfl.out.ch_demux_fastq

    EXTRACTUMI_wfl(ch_demux_fastq, ch_linkers)
    def ch_umi_fastq = EXTRACTUMI_wfl.out.ch_umi_fastq

    //changing umi output for FASTQCANDTRIM process. 
    def ch_umi_out = ch_umi_fastq
        .map { sample_id, fq1, fq2, log ->
            def meta = [ id: sample_id ]
            meta['single_end'] = (fq2 == null)
            tuple(meta, [fq1, fq2])
        }  
        // ch_ltr_out hould contain this type of data: [[id:D81_CART_d7_S3, single_end:false], [/beegfs/home/lrenteria/inspiired_nf/work/af/17177e890698137666321885a92ba0/D81_CART_d7_S3.ltr_filtered_R1.fastq.gz, /beegfs/home/lrenteria/inspiired_nf/work/af/17177e890698137666321885a92ba0/D81_CART_d7_S3.ltr_filtered_R2.fastq.gz]]
        //for debugging:  
        //      ch_ltr_out.view { "CH_LTR_OUT: ${it}" }

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


    LTRCHECKING_wfl(FASTQCANDTRIM_wfl.out.reads, ch_primer_ltr)
    def ch_ltr_chunks = LTRCHECKING_wfl.out.ch_ltr_chunks


    RCREMOVAL_wfl(ch_ltr_chunks, ch_linkers)
    //def ch_genomic_reads = RCREMOVAL_wfl.out.ch_rc_removed
        //DEBUGGING: we think that the reason why rcremoval is not caching is because: 
        //Nextflow is taking considerable time to cache the RCremoval outputs. 
        //When a subsequent process fails, the RCremoval outputs that have not finished caching are not counted as completed in the next execution round.
    def ch_all_genomic_reads = RCREMOVAL_wfl.out.ch_rc_removed.collect()
        // Split back into individual samples, to later send to findvectors
    def ch_split_samples = ch_all_genomic_reads.flatten()
        //DEBUGGING: 
        ch_split_samples.view{ "CH_SPLIT: ${ it } " }

    FINDVECTOR_wfl(ch_split_samples, ch_primer_ltr, file('vector.fasta'))
    def ch_vector_removed = FINDVECTOR_wfl.out.ch_vector_removed

    SHORTREMOVE_wfl(ch_vector_removed, ch_primer_ltr)
    def ch_short_removed = SHORTREMOVE_wfl.out.ch_short_removed

    DEREPLICATE_wfl(ch_short_removed)
    def ch_dereplicated = DEREPLICATE_wfl.out.ch_dereplicated

    // ALIGNMENT_wfl(ch_dereplicated, genome_name)

    //postprocessing workflow (processalignments)
            //in the postprocessing we will expand the reads with the keys file and divide the alignments in different groups. (unique, multihit, chimeric)
    //replicate/expand dereplicated reads

}