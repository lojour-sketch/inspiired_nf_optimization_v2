#!/usr/bin/env nextflow

// Parameters
params.samplesheet = '' // this sheet contains barcodes, processing parameters, LTR fragments, etc.
params.demuxSampleSheet = '' // this sheet contains only barcodes for demuxing
params.runfolderDir = ''
params.projectName = ''
params.FASTQfolderDir = ''


// Parameter validation
//if (!params.samplesheet) error "Missing --samplesheet parameter, provide path to samplesheet"
if (!params.runfolderDir) error "Missing --runfolderDir parameter. provide path to BCL data folder or tarball"
if (!params.projectName) error "Missing --projectName parameter, provide the name of the project to organize the result folders"
if (!params.FASTQfolderDir) error "Missing --FASTQfolderDir parameter, provide path to the folder with Undetermined FASTQ files"
if (!params.demuxSampleSheet) error "Missing --demuxSampleSheet parameter, provide path to the barcode samplesheet for demuxing"

// Log info

log.info "Using data from directory: ${params.runfolderDir}"

// Include modules

// MAIN

include { ALIGNMENT_wfl } from './subworkflows/local/alignment/main_insp_until_alignment'
include { POSTPROCESSING_twice_wfl } from './subworkflows/local/postprocessing/main_inspiired_insp_until_alignment'

// Workflow
workflow {

// ******************************** NECESSARY CHANNELS ********************************
        // we need a channel with the reference genome
        Channel.fromPath(params.samplesheet)
            .map { file ->
                def lines = file.text.readLines()
                def dataStart = lines.findIndexOf { it.startsWith('[Data]') } + 1
                def csvLines = lines[dataStart..-1].join("\n")
                return csvLines
            }
            .splitCsv(header: true)
            .map { row -> 
                def refGenomeMap = [
                    "hg19": "${params.runfolderDir}/../hg19_GRCh37_UCSC_initialrelease_2009.fa",
                    "hg38": "${params.runfolderDir}/../hg38_GRCh38_UCSC_initialrelease_2013.fa",
                    "hg18": "${params.runfolderDir}/../hg18_UCSC_2020_01_23_0222.fa"
                ]
                def refKnowngeneMap = [
                    "hg19": "TxDb.Hsapiens.UCSC.hg19.knownGene",
                    "hg38": "TxDb.Hsapiens.UCSC.hg38.refGene",
                    "hg18": "TxDb.Hsapiens.UCSC.hg18.knownGene"
                ]
                def refGenome = row.refGenome
                def refGenomeFile = refGenomeMap[refGenome]
                def refKnowngeneFile = refKnowngeneMap[refGenome]

                if (!refGenomeFile) {
                    error "Unsupported reference genome: ${refGenome}. Please use 'hg19' or 'hg38'."
                }
                return tuple(row.Sample_ID, refGenome, refGenomeFile, refKnowngeneFile)
            }
            .set { ch_refGenome }
        //we need a channel with the processing parameters
        Channel.fromPath(params.samplesheet)
            .map { file ->
                def lines = file.text.readLines()
                def dataStart = lines.findIndexOf { it.startsWith('[Data]') } + 1
                def csvLines = lines[dataStart..-1].join("\n")
                return csvLines
            }
            .splitCsv(header: true)
            .map { row ->
                tuple(row.Sample_ID, row.minPctIdent, row.maxAlignStart, row.maxFragLength)
            }
            .set { ch_processing_params }

// ************************************* WORKFLOW *************************************

    ALIGNMENT_wfl(ch_refGenome)
    def ch_aligned = ALIGNMENT_wfl.out.ch_aligned

    POSTPROCESSING_twice_wfl(ch_aligned, ch_processing_params, ch_refGenome)


}