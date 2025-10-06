//
// Read QC, UMI extraction and trimming
//

include { FASTQC_RAW           } from '../../../modules/nf-core/fastqc/main_raw'
include { UMITOOLS_EXTRACT } from '../../../modules/nf-core/umitools/extract/main'
include { TRIMGALORE       } from '../../../modules/nf-core/trimgalore/main'
include { FASTQC_TRIMMED    } from '../../../modules/nf-core/fastqc/main_trimmed'
// added a process to remove N containing sequences, in order to be able to do pairwise alignment later
include { N_REMOVING } from '../../../modules/local/nremoving_local/main'
//
// Function that parses TrimGalore log output file to get total number of reads after trimming
//
def getTrimGaloreReadsAfterFiltering(log_file) {
    def total_reads = 0
    def filtered_reads = 0
    log_file.eachLine { line ->
        def total_reads_matcher = line =~ /([\d\.]+)\ssequences processed in total/
        def filtered_reads_matcher = line =~ /shorter than the length cutoff[^:]+:\s*([\d\.]+)/
        if (total_reads_matcher) {
            total_reads = total_reads_matcher[0][1].toFloat()
        }
        if (filtered_reads_matcher) {
            filtered_reads = filtered_reads_matcher[0][1].toFloat()
        }
    }
    return total_reads - filtered_reads
}

workflow FASTQCANDTRIM_wfl {
    take:
    reads             // channel: [ val(meta), [ reads ] ]
    skip_fastqc       // boolean: true/false
    with_umi          // boolean: true/false
    skip_umi_extract  // boolean: true/false
    skip_trimming     // boolean: true/false
    umi_discard_read  // integer: 0, 1 or 2
    min_trimmed_reads // integer: > 0

    main:
    ch_versions = Channel.empty()
    fastqc_html = Channel.empty()
    fastqc_zip = Channel.empty()
    if (!skip_fastqc) {
        FASTQC_RAW(reads)
        fastqc_html = FASTQC_RAW.out.html
        fastqc_zip = FASTQC_RAW.out.zip
        ch_versions = ch_versions.mix(FASTQC_RAW.out.versions.first())
    }

    umi_reads = reads
    umi_log = Channel.empty()
    if (with_umi && !skip_umi_extract) {
        UMITOOLS_EXTRACT(reads)
        umi_reads = UMITOOLS_EXTRACT.out.reads
        umi_log = UMITOOLS_EXTRACT.out.log
        ch_versions = ch_versions.mix(UMITOOLS_EXTRACT.out.versions.first())

        // Discard R1 / R2 if required
        if (umi_discard_read in [1, 2]) {
            UMITOOLS_EXTRACT.out.reads
                .map { meta, reads ->
                    meta.single_end ? [meta, reads] : [meta + ['single_end': true], reads[umi_discard_read % 2]]
                }
                .set { umi_reads }
        }
    
    }

    trim_reads = umi_reads
    trim_unpaired = Channel.empty()
    trim_html = Channel.empty()
    trim_zip = Channel.empty()
    trim_log = Channel.empty()
    trim_read_count = Channel.empty()
    if (!skip_trimming) {
        TRIMGALORE(umi_reads)
        trim_unpaired = TRIMGALORE.out.unpaired
        trim_html = TRIMGALORE.out.html
        trim_zip = TRIMGALORE.out.zip
        trim_log = TRIMGALORE.out.log
        ch_versions = ch_versions.mix(TRIMGALORE.out.versions.first())

        //
        // Debugging
        //TRIMGALORE.out.trimmed_reads_N.view { meta, files -> 
        //    "Trimmed reads with N for sample ${meta}: ${files.join(', ')}"
        //}

       //TRIMGALORE.out.trimmed_reads_noN.view { meta, files -> 
        //    "Trimmed reads no N for sample ${meta}: ${files.join(', ')}"
        //}

        //
        // Filter FastQ files based on minimum trimmed read count after adapter trimming and N removing
        //

        TRIMGALORE.out.trimmed_reads
            .join(trim_log, remainder: true)
            .map { meta, reads_, trim_log_ ->
                if (trim_log) {
                    def num_reads = getTrimGaloreReadsAfterFiltering(meta.single_end ? trim_log_ : trim_log_[-1])
                    [meta, reads_, num_reads]
                }
                else {
                    [meta, reads, min_trimmed_reads.toFloat() + 1]
                }
            }
            .set { ch_num_trimmed_reads }

        ch_num_trimmed_reads
            .filter { _meta, _reads, num_reads -> num_reads >= min_trimmed_reads.toFloat() }
            .map { meta, reads_, _num_reads -> [meta, reads_] }
            .set { trim_reads }

        ch_num_trimmed_reads
            .map { meta, _reads, num_reads -> [meta, num_reads] }
            .set { trim_read_count }
            
    // Apart from filtering reads with trimgalore we also remove reads that contain Ns. we will later do a FASTQC on these reads
    N_REMOVING(trim_reads)
    trim_reads = N_REMOVING.out.reads

    fastqc_html_trimmed = Channel.empty()
    fastqc_zip_trimmed  = Channel.empty()
    if (!skip_trimming && !skip_fastqc) {
        FASTQC_TRIMMED(trim_reads)
        fastqc_html_trimmed = FASTQC_TRIMMED.out.html
        fastqc_zip_trimmed  = FASTQC_TRIMMED.out.zip
        ch_versions = ch_versions.mix(FASTQC_TRIMMED.out.versions.first())

    }
}

    emit:
    reads           = trim_reads // channel: [ val(meta), [ reads ] ]
    fastqc_html     // channel: [ val(meta), [ html ] ]
    fastqc_zip      // channel: [ val(meta), [ zip ] ]
    fastqc_html_trimmed     // channel: [ val(meta), [ html ] ]
    fastqc_zip_trimmed      // channel: [ val(meta), [ zip ] ]
    umi_log         // channel: [ val(meta), [ log ] ]
    trim_unpaired   // channel: [ val(meta), [ reads ] ]
    trim_html       // channel: [ val(meta), [ html ] ]
    trim_zip        // channel: [ val(meta), [ zip ] ]
    trim_log        // channel: [ val(meta), [ txt ] ]
    trim_read_count // channel: [ val(meta), val(count) ]
    versions        = ch_versions // channel: [ versions.yml ]
}

//for debugging
//      umi_reads must have this type of data: D83_CART_d0_S7, single_end=false, reads=[/beegfs/home/lrenteria/inspiired_nf/work/5e/664d31fe7f6f2e8a6303fe4b0ac925/D83_CART_d0_S7.umi_R1.fastq.gz, /beegfs/home/lrenteria/inspiired_nf/work/5e/664d31fe7f6f2e8a6303fe4b0ac925/D83_CART_d0_S7.umi_R2.fastq.gz]
//          umi_reads.view { "TRIMGALORE input tuple: ${it[0].id}, single_end=${it[0].single_end}, reads=${it[1]}" }