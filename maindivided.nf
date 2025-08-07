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
include { LTRREMOVING_wfl } from './subworkflows/local/LTR_removing/main'
include { RCREMOVAL_wfl } from './subworkflows/local/rcremoval/main'
include { SHORTREMOVE_wfl } from './subworkflows/local/short_seq_removal/main'

// Create necessary input tuple for bcl2fastq

Channel
    .of( tuple([id: 'run1'], file(params.samplesheet), file(params.runfolderDir, type: 'dir')) )
    .set { ch_bcl_input }



// Workflow
workflow {

    BCL2FASTQ_wfl(ch_bcl_input)
    def ch_demux_fastq = BCL2FASTQ_wfl.out.ch_demux_fastq

    EXTRACTUMI_wfl(ch_demux_fastq)
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

    LTRCHECKING_wfl(FASTQCANDTRIM_wfl.out.reads)
    def ch_ltr_checked = LTRCHECKING_wfl.out.ch_ltr_checked

    LTRREMOVING_wfl(ch_ltr_checked)
    def ch_ltr_removed = LTRREMOVING_wfl.out.ch_ltr_removed

    RCREMOVAL_wfl(ch_ltr_removed)
    def ch_rc_removed = RCREMOVAL_wfl.out.ch_rc_removed

    SHORTREMOVE_wfl(ch_rc_removed)
    def ch_short_removed = SHORTREMOVE_wfl.out.ch_short_removed

}