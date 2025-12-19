#!/usr/bin/env nextflow

// Parameters
params.samplesheet = '' // this sheet contains barcodes, processing parameters, LTR fragments, etc.
params.runfolderDir = ''
params.outputDir = ''
params.linkerdata = ''

// Parameter validation
//if (!params.samplesheet) error "Missing --samplesheet parameter, provide path to samplesheet"
if (!params.runfolderDir) error "Missing --runfolderDir p arameter. provide path to BCL data folder or tarball"
if (!params.outputDir) error "Missing --outputDir parameter, provide path to output directory"
if (!params.linkerdata) error "Missing --linkerdata parameter, provide .tsv file with linker sequence information"

// Log info

log.info "Using data from directory: ${params.runfolderDir}"

// Include modules

include { PREPROCESSING_wfl } from './subworkflows/local/preprocessing/blat_troubleshooting_main'
include { ALIGNMENT_wfl } from './subworkflows/local/alignment/blat_main'
//include { POSTPROCESSING_wfl } from './subworkflows/local/postprocessing/main'

// Create necessary input tuple for bcl2fastq. we only take two samples for faster troubleshooting

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
                tuple(row.Sample_ID, row.primer, row.ltrbit, row.largeLTRFrag, row.Sample_Project, row.mingDNA)
            }
            .set { ch_primer_ltr }
        // we need a channel with the reference genome
        Channel.fromPath(params.samplesheet)
            .map { file ->
                def lines = file.text.readLines()
                def dataStart = lines.findIndexOf { it.startsWith('[Data]') } + 1
                def csvLines = lines[dataStart..-1].join("\n")
                return csvLines
            }
            .splitCsv(header: true)
            .take(2)
            .map { row -> row.refGenome }
            .set { ch_refGenome }



// ************************************* WORKFLOW *************************************
    PREPROCESSING_wfl(ch_bcl_input, ch_linkers, ch_primer_ltr, file('vector.fasta') )
    def ch_dereplicated = PREPROCESSING_wfl.out.ch_dereplicated

    ALIGNMENT_wfl(ch_dereplicated, ch_refGenome)
   // def ch_aligned = ALIGNMENT_wfl.out.ch_aligned

    //POSTPROCESSING_wfl(ch_aligned, ch_primer_ltr)

            //in the postprocessing we will expand the reads with the keys file and divide the alignments in different groups. (unique, multihit, chimeric)
    //replicate/expand dereplicated reads

}