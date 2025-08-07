include { BCL2FASTQ_local } from '/home/lrenteria/inspiired_nf/modules/local/bcl2fastq_local/main.nf'

workflow BCL2FASTQ_wfl {
    take:
    ch_bcl_input


    main:
    // Run bcl2fastq process
    BCL2FASTQ_local(ch_bcl_input)

    //save bcl2fastq output channels

        BCL2FASTQ_local.out.fastq           .collect().set { ch_demux_fastq }
        BCL2FASTQ_local.out.fastq_idx       .collect().set { ch_fastq_idx }
        BCL2FASTQ_local.out.undetermined    .collect().set { ch_undetermined }
        BCL2FASTQ_local.out.undetermined_idx.collect().set { ch_undetermined_idx }
        BCL2FASTQ_local.out.reports         .collect().set { ch_bcl_reports }
        BCL2FASTQ_local.out.stats           .collect().set { ch_bcl_stats }
        BCL2FASTQ_local.out.interop         .collect().set { ch_bcl_interop }
        BCL2FASTQ_local.out.versions        .collect().set { ch_bcl_versions }


    emit:
    ch_demux_fastq
    ch_fastq_idx
    ch_undetermined
    ch_undetermined_idx
    ch_bcl_reports
    ch_bcl_stats
    ch_bcl_interop
    ch_bcl_versions
}

//for debugging:
    //ch_demux_fastq must have this type of data: [[id:run1], [/beegfs/home/lrenteria/inspiired_nf/work/5c/2e29e2c2b7ed467ace2fa2d81ffe71/results/alloCART/D81_CART-KO_d07_S4_R1_001.fastq.gz, /beegfs/home/lrenteria/inspiired_nf/work/5c/2e29e2c2b7ed467ace2fa2d81ffe71/results/alloCART/D81_CART-KO_d07_S4_R2_001.fastq.gz, /beegfs/home/lrenteria/inspiired_nf/work/5c/2e29e2c2b7ed467ace2fa2d81ffe71/results/alloCART/D81_CART-KO_d0_S2_R1_001.fastq.gz, /beegfs/home/lrenteria/inspiired_nf/work/5c/2e29e2c2b7ed467ace2fa2d81ffe71/results/alloCART/D81_CART-KO_d0_S2_R2_001.fastq.gz, /beegfs/home/lrenteria/inspiired_nf/work/5c/2e29e2c2b7ed467ace2fa2d81ffe71/results/alloCART/D81_CART-KO_d14_S6_R1_001.fastq.gz, /beegfs/home/lrenteria/inspiired_nf/work/5c/2e29e2c2b7ed467ace2fa2d81ffe71/results/alloCART/D81_CART-KO_d14_S6_R2_001.fastq.gz, /beegfs/home/lrenteria/inspiired_nf/work/5c/2e29e2c2b7ed467ace2fa2d81ffe71/results/alloCART/D81_CART_d0_S1_R1_001.fastq.gz, /beegfs/home/lrenteria/inspiired_nf/work/5c/2e29e2c2b7ed467ace2fa2d81ffe71/results/alloCART/D81_CART_d0_S1_R2_001.fastq.gz, /beegfs/home/lrenteria/inspiired_nf/work/5c/2e29e2c2b7ed467ace2fa2d81ffe71/results/alloCART/D81_CART_d14_S5_R1_001.fastq.gz, /beegfs/home/lrenteria/inspiired_nf/work/5c/2e29e2c2b7ed467ace2fa2d81ffe71/results/alloCART/D81_CART_d14_S5_R2_001.fastq.gz, /beegfs/home/lrenteria/inspiired_nf/work/5c/2e29e2c2b7ed467ace2fa2d81ffe71/results/alloCART/D81_CART_d7_S3_R1_001.fastq.gz, /beegfs/home/lrenteria/inspiired_nf/work/5c/2e29e2c2b7ed467ace2fa2d81ffe71/results/alloCART/D81_CART_d7_S3_R2_001.fastq.gz, /beegfs/home/lrenteria/inspiired_nf/work/5c/2e29e2c2b7ed467ace2fa2d81ffe71/results/alloCART/D83_CART-KO_d07_S10_R1_001.fastq.gz, /beegfs/home/lrenteria/inspiired_nf/work/5c/2e29e2c2b7ed467ace2fa2d81ffe71/results/alloCART/D83_CART-KO_d07_S10_R2_001.fastq.gz, /beegfs/home/lrenteria/inspiired_nf/work/5c/2e29e2c2b7ed467ace2fa2d81ffe71/results/alloCART/D83_CART-KO_d0_S8_R1_001.fastq.gz, /beegfs/home/lrenteria/inspiired_nf/work/5c/2e29e2c2b7ed467ace2fa2d81ffe71/results/alloCART/D83_CART-KO_d0_S8_R2_001.fastq.gz, /beegfs/home/lrenteria/inspiired_nf/work/5c/2e29e2c2b7ed467ace2fa2d81ffe71/results/alloCART/D83_CART-KO_d14_S12_R1_001.fastq.gz, /beegfs/home/lrenteria/inspiired_nf/work/5c/2e29e2c2b7ed467ace2fa2d81ffe71/results/alloCART/D83_CART-KO_d14_S12_R2_001.fastq.gz, /beegfs/home/lrenteria/inspiired_nf/work/5c/2e29e2c2b7ed467ace2fa2d81ffe71/results/alloCART/D83_CART_d0_S7_R1_001.fastq.gz, /beegfs/home/lrenteria/inspiired_nf/work/5c/2e29e2c2b7ed467ace2fa2d81ffe71/results/alloCART/D83_CART_d0_S7_R2_001.fastq.gz, /beegfs/home/lrenteria/inspiired_nf/work/5c/2e29e2c2b7ed467ace2fa2d81ffe71/results/alloCART/D83_CART_d14_S11_R1_001.fastq.gz, /beegfs/home/lrenteria/inspiired_nf/work/5c/2e29e2c2b7ed467ace2fa2d81ffe71/results/alloCART/D83_CART_d14_S11_R2_001.fastq.gz, /beegfs/home/lrenteria/inspiired_nf/work/5c/2e29e2c2b7ed467ace2fa2d81ffe71/results/alloCART/D83_CART_d7_S9_R1_001.fastq.gz, /beegfs/home/lrenteria/inspiired_nf/work/5c/2e29e2c2b7ed467ace2fa2d81ffe71/results/alloCART/D83_CART_d7_S9_R2_001.fastq.gz]]
    //        BCL2FASTQ_local.out.fastq.collect().view { "Demultiplexed FASTQ files: ${it}" }.set { ch_demux_fastq }
    //ch_fastq_idx must have this type of data: [[id:run1], [/beegfs/home/lrenteria/inspiired_nf/work/5c/2e29e2c2b7ed467ace2fa2d81ffe71/results/alloCART/D81_CART-KO_d07_S4_I1_001.fastq.gz, /beegfs/home/lrenteria/inspiired_nf/work/5c/2e29e2c2b7ed467ace2fa2d81ffe71/results/alloCART/D81_CART-KO_d07_S4_I2_001.fastq.gz, /beegfs/home/lrenteria/inspiired_nf/work/5c/2e29e2c2b7ed467ace2fa2d81ffe71/results/alloCART/D81_CART-KO_d0_S2_I1_001.fastq.gz, /beegfs/home/lrenteria/inspiired_nf/work/5c/2e29e2c2b7ed467ace2fa2d81ffe71/results/alloCART/D81_CART-KO_d0_S2_I2_001.fastq.gz, /beegfs/home/lrenteria/inspiired_nf/work/5c/2e29e2c2b7ed467ace2fa2d81ffe71/results/alloCART/D81_CART-KO_d14_S6_I1_001.fastq.gz, /beegfs/home/lrenteria/inspiired_nf/work/5c/2e29e2c2b7ed467ace2fa2d81ffe71/results/alloCART/D81_CART-KO_d14_S6_I2_001.fastq.gz, /beegfs/home/lrenteria/inspiired_nf/work/5c/2e29e2c2b7ed467ace2fa2d81ffe71/results/alloCART/D81_CART_d0_S1_I1_001.fastq.gz, /beegfs/home/lrenteria/inspiired_nf/work/5c/2e29e2c2b7ed467ace2fa2d81ffe71/results/alloCART/D81_CART_d0_S1_I2_001.fastq.gz, /beegfs/home/lrenteria/inspiired_nf/work/5c/2e29e2c2b7ed467ace2fa2d81ffe71/results/alloCART/D81_CART_d14_S5_I1_001.fastq.gz, /beegfs/home/lrenteria/inspiired_nf/work/5c/2e29e2c2b7ed467ace2fa2d81ffe71/results/alloCART/D81_CART_d14_S5_I2_001.fastq.gz, /beegfs/home/lrenteria/inspiired_nf/work/5c/2e29e2c2b7ed467ace2fa2d81ffe71/results/alloCART/D81_CART_d7_S3_I1_001.fastq.gz, /beegfs/home/lrenteria/inspiired_nf/work/5c/2e29e2c2b7ed467ace2fa2d81ffe71/results/alloCART/D81_CART_d7_S3_I2_001.fastq.gz, /beegfs/home/lrenteria/inspiired_nf/work/5c/2e29e2c2b7ed467ace2fa2d81ffe71/results/alloCART/D83_CART-KO_d07_S10_I1_001.fastq.gz, /beegfs/home/lrenteria/inspiired_nf/work/5c/2e29e2c2b7ed467ace2fa2d81ffe71/results/alloCART/D83_CART-KO_d07_S10_I2_001.fastq.gz, /beegfs/home/lrenteria/inspiired_nf/work/5c/2e29e2c2b7ed467ace2fa2d81ffe71/results/alloCART/D83_CART-KO_d0_S8_I1_001.fastq.gz, /beegfs/home/lrenteria/inspiired_nf/work/5c/2e29e2c2b7ed467ace2fa2d81ffe71/results/alloCART/D83_CART-KO_d0_S8_I2_001.fastq.gz, /beegfs/home/lrenteria/inspiired_nf/work/5c/2e29e2c2b7ed467ace2fa2d81ffe71/results/alloCART/D83_CART-KO_d14_S12_I1_001.fastq.gz, /beegfs/home/lrenteria/inspiired_nf/work/5c/2e29e2c2b7ed467ace2fa2d81ffe71/results/alloCART/D83_CART-KO_d14_S12_I2_001.fastq.gz, /beegfs/home/lrenteria/inspiired_nf/work/5c/2e29e2c2b7ed467ace2fa2d81ffe71/results/alloCART/D83_CART_d0_S7_I1_001.fastq.gz, /beegfs/home/lrenteria/inspiired_nf/work/5c/2e29e2c2b7ed467ace2fa2d81ffe71/results/alloCART/D83_CART_d0_S7_I2_001.fastq.gz, /beegfs/home/lrenteria/inspiired_nf/work/5c/2e29e2c2b7ed467ace2fa2d81ffe71/results/alloCART/D83_CART_d14_S11_I1_001.fastq.gz, /beegfs/home/lrenteria/inspiired_nf/work/5c/2e29e2c2b7ed467ace2fa2d81ffe71/results/alloCART/D83_CART_d14_S11_I2_001.fastq.gz, /beegfs/home/lrenteria/inspiired_nf/work/5c/2e29e2c2b7ed467ace2fa2d81ffe71/results/alloCART/D83_CART_d7_S9_I1_001.fastq.gz, /beegfs/home/lrenteria/inspiired_nf/work/5c/2e29e2c2b7ed467ace2fa2d81ffe71/results/alloCART/D83_CART_d7_S9_I2_001.fastq.gz]]
    //        BCL2FASTQ_local.out.fastq_idx.collect().view { "Demultiplexed FASTQ index files: ${it}" }.set { ch_fastq_idx }
    //ch_undetermuned must have this type of data: [[id:run1], [/beegfs/home/lrenteria/inspiired_nf/work/5c/2e29e2c2b7ed467ace2fa2d81ffe71/results/Undetermined_S0_R1_001.fastq.gz, /beegfs/home/lrenteria/inspiired_nf/work/5c/2e29e2c2b7ed467ace2fa2d81ffe71/results/Undetermined_S0_R2_001.fastq.gz]]
    //        BCL2FASTQ_local.out.undetermined.collect().view { "Undetermined FASTQ files: ${it}" }.set { ch_undetermined }
    //ch_undetermined_idx must have this type of data: [[id:run1], [/beegfs/home/lrenteria/inspiired_nf/work/5c/2e29e2c2b7ed467ace2fa2d81ffe71/results/Undetermined_S0_I1_001.fastq.gz, /beegfs/home/lrenteria/inspiired_nf/work/5c/2e29e2c2b7ed467ace2fa2d81ffe71/results/Undetermined_S0_I2_001.fastq.gz]]
    //        BCL2FASTQ_local.out.undetermined_idx.collect().view { "Undetermined FASTQ index files: ${it}" }.set { ch_undetermined_idx }
    //For other created channels:
    //        BCL2FASTQ_local.out.reports         .collect().view { "BCL2FASTQ reports: ${it}" }.set { ch_bcl_reports }
    //        BCL2FASTQ_local.out.stats           .collect().view { "BCL2FASTQ stats: ${it}" }.set { ch_bcl_stats }
    //        BCL2FASTQ_local.out.interop         .collect().view { "BCL2FASTQ interop files: ${it}" }.set { ch_bcl_interop }
    //        BCL2FASTQ_local.out.versions        .collect().view { "BCL2FASTQ versions: ${it}" }.set { ch_bcl_versions }