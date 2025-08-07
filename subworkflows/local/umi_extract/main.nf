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

        // ch_reads_by_sample must have this type of data: D83_CART_d14_S11
        //for debugging:
        //ch_reads_by_sample.view { "CH_READS_BY_SAMPLE: ${it}" }

    //create channel containing linker sequeces that constrain the UMI sequence
    ch_linkers = Channel
    .fromPath(params.linkerdata)
    .splitCsv(header: true)
    .map { row ->
        tuple(row.sample_id, row.sample_unique_linker, row.common_linker)
    }

        // ch_linkers must have this type of data: [D81_CART-KO_d07_S4, CGGCTTACAATTCCTGCGAC, CTCCGCTTAAGGGACT]
        //for debugging:
        //ch_linkers.view { "CH_LINKERS channel: ${it}" }
       

    //now we combine the reads channel with the linkers channel per sample
    ch_reads_by_sample
        .join(ch_linkers) //join by sample_id
        .map { sample_id, reads, linker1, linker2 ->
                tuple(sample_id, linker1, linker2, reads)
        }
        .set { ch_umi_extract_input }

        //ch_umi_extract_input must have this type of data: [D83_CART_d7_S9, GAACGAGCACTAGTAAGCCC, CTCCGCTTAAGGGACT, [/beegfs/home/lrenteria/inspiired_nf/work/5c/2e29e2c2b7ed467ace2fa2d81ffe71/results/alloCART/D83_CART_d7_S9_R1_001.fastq.gz, /beegfs/home/lrenteria/inspiired_nf/work/5c/2e29e2c2b7ed467ace2fa2d81ffe71/results/alloCART/D83_CART_d7_S9_R2_001.fastq.gz]]
        //for debugging:
        //ch_umi_extract_input.view{ "CH_UMI_EXTRACT_INPUT: ${it}" }

    //now we run the umi extract in only R1 files
    UMI_EXTRACT_LOCAL(ch_umi_extract_input)

    //tupleoutput channel contains: 
    UMI_EXTRACT_LOCAL.out
    .set { ch_umi_fastq }

        //ch_umi_fastq must have this type of data: [D81_CART_d7_S3, /beegfs/home/lrenteria/inspiired_nf/work/bc/c836f1144ba7d6192ad2e569dacecf/D81_CART_d7_S3.umi_R1.fastq.gz, /beegfs/home/lrenteria/inspiired_nf/work/bc/c836f1144ba7d6192ad2e569dacecf/D81_CART_d7_S3.umi_R2.fastq.gz, /beegfs/home/lrenteria/inspiired_nf/work/bc/c836f1144ba7d6192ad2e569dacecf/D81_CART_d7_S3.umi_extract.log]
        //for debugging:
        //ch_umi_fastq.view{ "CH_UMI_FASTQ: ${it}" }

    emit:
    ch_umi_fastq

}