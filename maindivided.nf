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
include { EXTRACTUMI_wfl } from './subworkflows/local/umi_extract_local/main.nf'
include { FASTQCANDTRIM_wfl } from './subworkflows/nf-core/fastq_fastqc_umitools_trimgalore/main.nf'
//include { MULTIQC } from './subworkflows/nf-core/multiqc/main.nf'

// Create necessary input tuple for bcl2fastq

Channel
    .of( tuple([id: 'run1'], file(params.samplesheet), file(params.runfolderDir, type: 'dir')) )
    .set { ch_bcl_input }



// Workflow
workflow {

    BCL2FASTQ_wfl(ch_bcl_input)
    EXTRACTUMI_wfl(BCL2FASTQ_wfl.out.ch_demux_fastq)

    //changing umi output for FASTQCANDTRIM process
    ch_umi_out = EXTRACTUMI_wfl.out.ch_umi_fastq
        .map { sample_id, fq1, fq2, log ->
            def meta = [ id: sample_id ]
            meta['single_end'] = (fq2 == null)
            tuple(meta, [fq1, fq2])
        }
        .view{"FASTQC input tuple: ${it}"}    

    FASTQCANDTRIM_wfl(
        ch_umi_out,           // input tuple: [ meta, reads ]
        false,                  //skip_fastqc
        false,                  // with_umi
        true,                   //skip_umi_extract
        false,                  // skip_trimming
        0,                   // umi_discard_read (0 , 1 or 2)
        10000,                   // min_trimmed_reads (default 10000)
    )


}