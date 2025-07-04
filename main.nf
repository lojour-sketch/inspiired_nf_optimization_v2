#!/usr/bin/env nextflow

// Parameters
params.samplesheet = ''
params.runfolderDir = ''
params.outputDir = ''

// Parameter validation
if (!params.samplesheet) error "Missing --samplesheet parameter, provide path to samplesheet"
if (!params.runfolderDir) error "Missing --runfolderDir parameter. provide path to BCL data folder or tarball"
if (!params.outputDir) error "Missing --outputDir parameter, provide path to output directory"

// Log info

log.info "Using data from directory: ${params.runfolderDir}"

// Include modules

include { BCL2FASTQ } from './modules/nf-core/bcl2fastq/main.nf' 
include { BCL2FASTQ_local } from './modules/local/bcl2fastq_local/main.nf'

// Create necessary input tuple for bcl2fastq

Channel
    .of( tuple([id: 'run1'], file(params.samplesheet), file(params.runfolderDir)) )
    .set { ch_bcl_input }

// Workflow
workflow {
    log.info "\n"
    log.info "************* Starting BCL conversion to FASTQ ****************"
    log.info "\n"

    // Run bcl2fastq process
    BCL2FASTQ_local(ch_bcl_input)

    //save bcl2fastq output channels

    BCL2FASTQ_local.out.fastq           .set { ch_demux_fastq }
    BCL2FASTQ_local.out.fastq_idx       .set { ch_fastq_idx }
    BCL2FASTQ_local.out.undetermined    .set { ch_undetermined }
    BCL2FASTQ_local.out.undetermined_idx.set { ch_undetermined_idx }
    BCL2FASTQ_local.out.reports         .set { ch_reports }
    BCL2FASTQ_local.out.stats           .set { ch_stats }
    BCL2FASTQ_local.out.interop         .set { ch_interop }
    BCL2FASTQ_local.out.versions        .set { ch_versions }


}