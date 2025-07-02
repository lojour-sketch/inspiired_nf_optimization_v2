#!/usr/bin/env nextflow

// Parameters
params.samplesheet = 'samplesheet.csv'
params.bcldata = 'bcl_data.tar.gz'
params.output = 'results'

// Parameter validation
if (!params.samplesheet) error "Missing --samplesheet parameter, provide path to samplesheet"
if (!params.bcldata) error "Missing --bcldata parameter. provide path to BCL data folder or tarball"
if (!params.output) error "Missing --output parameter, provide path to output directory"

// Log info

log.info "Starting pipeline for job ID: ${params.jobid}"
log.info "Using data from directory: ${params.input}"

// Include modules

include { BCL2FASTQ } from '../modules/nf-core/bcl2fastq/main' 

// Create necessary input tuple for bcl2fastq

Channel
    .of( tuple([id: 'run1'], file(params.samplesheet), file(params.bcldata)) )
    .set { ch_bcl_input }

// Workflow
workflow {
    log.info "\n"
    log.info "************* Starting BCL conversion to FASTQ ****************"
    log.info "\n"

    BCLFASTQ(ch_bcl_input)
}