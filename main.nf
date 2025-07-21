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

//include { BCL2FASTQ } from './modules/nf-core/bcl2fastq/main.nf' 
include { BCL2FASTQ_local } from './modules/local/bcl2fastq_local/main.nf'
include { UMI_EXTRACT_LOCAL } from './modules/local/umi_extract_local/main.nf'
include { FASTQ_FASTQC_UMITOOLS_TRIMGALORE } from './subworkflows/nf-core/fastq_fastqc_umitools_trimgalore/main.nf'
include { MULTIQC } from './modules/nf-core/multiqc/main.nf'

// Create necessary input tuple for bcl2fastq

Channel
    .of( tuple([id: 'run1', single_end: false], file(params.samplesheet), file(params.runfolderDir, type: 'dir')) )
    .set { ch_bcl_input }



// Workflow
workflow {
    log.info "\n"
    log.info "************* Starting BCL conversion to FASTQ ****************"
    log.info "\n"

    // Run bcl2fastq process
    BCL2FASTQ_local(ch_bcl_input)

    //save bcl2fastq output channels

    BCL2FASTQ_local.out.fastq           .collect().view { "Demultiplexed FASTQ files: ${it}" }.set { ch_demux_fastq }
    BCL2FASTQ_local.out.fastq_idx       .collect().view { "Demultiplexed FASTQ index files: ${it}" }.set { ch_fastq_idx }
    BCL2FASTQ_local.out.undetermined    .collect().view { "Undetermined FASTQ files: ${it}" }.set { ch_undetermined }
    BCL2FASTQ_local.out.undetermined_idx.collect().view { "Undetermined FASTQ index files: ${it}" }.set { ch_undetermined_idx }
    BCL2FASTQ_local.out.reports         .collect().view { "BCL2FASTQ reports: ${it}" }.set { ch_reports }
    BCL2FASTQ_local.out.stats           .collect().view { "BCL2FASTQ stats: ${it}" }.set { ch_stats }
    BCL2FASTQ_local.out.interop         .collect().view { "BCL2FASTQ interop files: ${it}" }.set { ch_interop }
    BCL2FASTQ_local.out.versions        .collect().view { "BCL2FASTQ versions: ${it}" }.set { ch_versions }


    //crete necessary input channel for umi_extract_local
    //first create channel with tuple(sample_id, [R1, R2])

    ch_reads_by_sample = ch_demux_fastq
        .map { id, file_list -> file_list }  // ignore id (e.g., "run1")
        .flatten()
        .map { file ->
            def name = file.getFileName().toString()
            def sample_name = name.replaceAll(/_R[12]_001\.fastq\.gz$/, '')
            tuple(sample_name, file)
        }
        .groupTuple()
        .view { sample_id, pair -> 
            "FASTQ pair: ${sample_id} -> ${pair}"
        }

    //create channel containing linker sequeces that constrain the UMI sequence
    ch_linkers = Channel
    .fromPath(params.linkerdata)
    .splitCsv(header: true)
    .map { row ->
        tuple(row.sample_id, row.sample_unique_linker, row.common_linker)
    }

    ch_linkers.view { println "ch_linkers: $it" }


    //debugging
    ch_reads_by_sample.view { "READS SAMPLES: ${it[0]}" }
    ch_linkers
    .view { id, l1, l2 -> "LINKER tuple: ${id} -> ${l1}, ${l2}" }
    
    ch_linkers.view { "LINKERS SAMPLES: ${it[0]}" }

    ch_linkers.view { "Linker full input: ${it}" }


    //now we combine the reads channel with the linkers channel per sample
    ch_reads_by_sample
        .join(ch_linkers) //join by sample_id
        .map { sample_id, reads, linker1, linker2 ->
                tuple(sample_id, linker1, linker2, reads)
        }
        .view { "Reads and linkers joined channel: ${it}"}
        .set { ch_umi_extract_input }

    //now we run the umi extract in only R1 files
    log.info "\n"
    log.info "************* Starting UMI extraction ****************"
    log.info "\n"
    UMI_EXTRACT_LOCAL(ch_umi_extract_input)

    UMI_EXTRACT_LOCAL.out
    .view { "UMI extraction output: ${it}" }
    .set { ch_umi_fastq }


    //change umi output to match the next process' input: tuple: [ meta, reads ]
    ch_umi_out = ch_umi_fastq
        .map { sample_id, fq1, fq2, log ->
            def meta = [ id: sample_id ]
            meta['single_end'] = (fq2 == null)
            tuple(meta, [fq1, fq2])
        }
        .view{"FASTQC input tuple: ${it}"}

    //perform FASTQ_FASTQC_UMITOOLS_TRIMGALORE in each sample's reads 
    //we check the bases quality, extract umi sequences and trim low quality bases

    log.info "\n"
    log.info "************* Starting FASTQC analysis and trimming ****************"
    log.info "\n"
    FASTQ_FASTQC_UMITOOLS_TRIMGALORE(
        ch_umi_out,           // input tuple: [ meta, reads ]
        false,                  //skip_fastqc
        false,                  // with_umi
        true,                   //skip_umi_extract
        false,                  // skip_trimming
        0,                   // umi_discard_read ()
        10000,                   // min_trimmed_reads (default 10000)
    )

    FASTQ_FASTQC_UMITOOLS_TRIMGALORE.out.reads
        .view { "Trimmed paired reads: ${it}" }
        .set { ch_trimmed_reads }

    FASTQ_FASTQC_UMITOOLS_TRIMGALORE.out.fastqc_html
        .view { "FASTQC HTML reports (raw): ${it}" }
        .set { ch_fastqc_html }

    FASTQ_FASTQC_UMITOOLS_TRIMGALORE.out.fastqc_zip
        .view { "FASTQC ZIP reports (raw): ${it}" }
        .set { ch_fastqc_zip }

    FASTQ_FASTQC_UMITOOLS_TRIMGALORE.out.umi_log
        .view { "UMI-tools extract logs: ${it}" }
        .set { ch_umi_log }

    FASTQ_FASTQC_UMITOOLS_TRIMGALORE.out.trim_unpaired
        .view { "Unpaired reads after trimming: ${it}" }
        .set { ch_trim_unpaired }

    FASTQ_FASTQC_UMITOOLS_TRIMGALORE.out.trim_html
        .view { "TrimGalore FASTQC HTML: ${it}" }
        .set { ch_trim_html }

    FASTQ_FASTQC_UMITOOLS_TRIMGALORE.out.trim_zip
        .view { "TrimGalore FASTQC ZIP: ${it}" }
        .set { ch_trim_zip }

    FASTQ_FASTQC_UMITOOLS_TRIMGALORE.out.trim_log
        .view { "TrimGalore logs: ${it}" }
        .set { ch_trim_log }

    FASTQ_FASTQC_UMITOOLS_TRIMGALORE.out.trim_read_count
        .view { "Trimmed read count: ${it}" }
        .set { ch_trim_read_count }

    FASTQ_FASTQC_UMITOOLS_TRIMGALORE.out.versions
        .view { "Software versions: ${it}" }
        .set { ch_versions }

    //MULTIQC(FASTQ_FASTQC_UMITOOLS_TRIMGALORE.out.fastqc_zip)
    //MULTIQC(FASTQ_FASTQC_UMITOOLS_TRIMGALORE.out.trim_zip)

}