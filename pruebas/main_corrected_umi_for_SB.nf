#!/usr/bin/env nextflow

// Parameters
params.BCLorFASTQ = '' //options: BCL or FASTQ
params.runfolderDir = ''
params.projectName = ''
params.samplesheet = '' // this sheet contains barcodes, processing parameters, LTR fragments, etc.

params.FASTQfolderDir = '' // path to the folder with Undetermined FASTQ files (only for FASTQ input)
params.readStructure = '' // this is the read structure of the FASTQ files (only for FASTQ input)
params.demuxSampleSheet = '' // this sheet contains only barcodes for demuxing (only for FASTQ input)



// Parameter validation
//if (!params.samplesheet) error "Missing --samplesheet parameter, provide path to samplesheet"
if (!params.runfolderDir) error "Missing --runfolderDir parameter. provide path to BCL data folder or tarball"
if (!params.projectName) error "Missing --projectName parameter, provide the name of the project to organize the result folders"
if (!params.BCLorFASTQ) error "Missing --BCLorFASTQ parameter, provide the type of data: BCL or FASTQ"
if (params.BCLorFASTQ != 'BCL' && params.BCLorFASTQ != 'FASTQ')  error "--BCLorFASTQ parameter must be either BCL or FASTQ"
if (params.BCLorFASTQ == 'FASTQ') {
    if (!params.FASTQfolderDir) error "Missing --FASTQfolderDir parameter, provide path to the folder with Undetermined FASTQ files"
    if (!params.readStructure) error "Missing --readStructure parameter, provide the read structure of the FASTQ files"
    if (!params.demuxSampleSheet) error "Missing --demuxSampleSheet parameter, provide path to the barcode samplesheet for demuxing"
}
if (!params.samplesheet) error "Missing --samplesheet parameter, provide path to samplesheet"



// Log info

log.info "Using data from directory: ${params.runfolderDir}"

// Include modules

// MAIN

if (params.BCLorFASTQ == 'BCL') {
    include { PREPROCESSING_wfl } from './subworkflows/local/preprocessing/main_bcl'
} else {
    include { PREPROCESSING_wfl } from './subworkflows/local/preprocessing/main_fastq_SB'
}
include { ALIGNMENT_wfl } from './subworkflows/local/alignment/main'
include { POSTPROCESSING_wfl } from './subworkflows/local/postprocessing/main'


if (params.BCLorFASTQ == 'BCL') {
    // Create necessary input tuple for bcl2fastq

    Channel
        .of( tuple([id: 'run1'], file(params.samplesheet), file(params.runfolderDir, type: 'dir')) )
        .set { ch_bcl_input }
}

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
                return tuple(refGenome, refGenomeFile, refKnowngeneFile)
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
                tuple(row.minPctIdent, row.maxAlignStart, row.maxFragLength)
            }
            .set { ch_processing_params }
        // we need a channel with the vector sequence
        Channel.fromPath(params.samplesheet)
            .map { file ->
                def lines = file.text.readLines()
                def dataStart = lines.findIndexOf { it.startsWith('[Data]') } + 1
                def csvLines = lines[dataStart..-1].join("\n")
                return csvLines
            }
            .splitCsv(header: true)
            .map { row ->
                def vectorSeq = row.vectorSeq
                return tuple("${params.runfolderDir}/${vectorSeq}")
            }
            .set { ch_vector }

// ************************************* WORKFLOW *************************************
    if (params.BCLorFASTQ == 'BCL') {
        PREPROCESSING_wfl(ch_bcl_input, ch_linkers, ch_primer_ltr, ch_vector )
    } else {
        PREPROCESSING_wfl(ch_linkers, ch_primer_ltr, ch_vector )
    }
    def ch_short_removed = PREPROCESSING_wfl.out.ch_short_removed

    ALIGNMENT_wfl(ch_short_removed, ch_refGenome)
    def ch_aligned = ALIGNMENT_wfl.out.ch_aligned

    POSTPROCESSING_wfl(ch_aligned, ch_processing_params, ch_refGenome)


}