include { UMI_EXTRACT_LOCAL } from '/home/lrenteria/inspiired_nf/modules/local/umi_extract_local/main.nf'

workflow EXTRACTUMI_wfl {
    take:
    ch_demux_fastq

    main:
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

    emit:
    ch_umi_fastq

}