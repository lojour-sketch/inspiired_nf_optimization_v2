include { LTRchecking_seqkit } from '../../../modules/local/LTR_primer_checking_local/main_seqkit_fullfq'
//include { splitFastq } from '../../../modules/local/splitmergefastq/splitmain'
//include { mergeChunks } from '../../../modules/local/splitmergefastq/mergemain'

workflow LTRCHECKING_wfl {
    take:
    reads
    primerltr

    main:

    reads
        .map { meta, reads -> tuple(meta.id, reads) }
        .join( primerltr.map { id, primer, ltrbit, largeLTRFrag, mingDNA -> tuple(id, primer, ltrbit, largeLTRFrag, mingDNA) } )
        .map { sample_id, reads, primer, ltrbit, largeLTRFrag, mingDNA ->
            tuple(sample_id, reads, primer, ltrbit, largeLTRFrag, mingDNA)
        }
        .set { ch_input }

    //following lines for when we want to do it with splitted fastqs
    // splitFastq(ch_input)
    // splitFastq.out.view { "Split fastq output ${it}" }


    // //prepare input for LTRchecking_seqkit
    // //we are assuming that chunksr1 and chunksr2 have matching files by name and order. 
    // //we need to do this to give path(reads) (but in different chunks to save memory) to the LTRchecking_seqkit process

    // splitFastq.out
    //     .map { meta, chunks_r1, chunks_r2, primer, ltrbit, largeLTRFrag, mingDNA ->
    //         def sample_id = meta instanceof Map ? meta.id : meta
    //         tuple(sample_id, chunks_r1, chunks_r2, primer, ltrbit, largeLTRFrag, mingDNA)
    //     }
    //     .flatMap { sample_id, chunks_r1, chunks_r2, primer, ltrbit, largeLTRFrag, mingDNA ->
    //         def chunk_pairs = []
    //         def r1_files = chunks_r1.toList().sort()
    //         def r2_files = chunks_r2.toList().sort()
    //         r1_files.eachWithIndex { f1, i ->
    //             def f2 = r2_files[i]
    //             if (!f2) throw new IllegalStateException("No matching R2 chunk for ${f1}")
    //             def matcher = (f1.name =~ /part_(\d+)\.fq\.gz/)
    //             def idx = matcher ? matcher[0][1] : i
    //             chunk_pairs << tuple(sample_id, f1, f2, primer, ltrbit, largeLTRFrag, mingDNA, idx)
    //         }
    //         return chunk_pairs
    //     }
    //     .view { "LTR checking input: ${it}" }
    //     .set { ch_LTRchecking_seqkit_input }

    LTRchecking_seqkit(ch_input)
    // LTRchecking_seqkit.out.view { "LTR checking output: ${it}" }.set{ ch_ltr_chunks }
    LTRchecking_seqkit.out.reads.set{ ch_ltr_chunks } //[D83_CART_d7_S9, /beegfs/home/lrenteria/inspiired_nf/work/a4/6c2ff6c57d7ef130a7b82a3e4b27eb/D83_CART_d7_S9.ltr_filtered_R1.chunk004.fastq.gz, /beegfs/home/lrenteria/inspiired_nf/work/a4/6c2ff6c57d7ef130a7b82a3e4b27eb/D83_CART_d7_S9.ltr_filtered_R2.chunk004.fastq.gz, GAAATC, TCTAGCA, CAGACCCTTTTAGTCAGTGTGGAAAATCTCTAGCA, 30, 004]

    

    emit:
    ch_ltr_chunks      //the format is: [sample, LTRfilteredR1, LTRfilteredR2, primer, ltrbit, largeLTRFrag, mingDNA, idx]

}