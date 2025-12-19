include { INDEX_SORT_BAM_local } from '../../../modules/local/index_sort_bam_local/main'
include { BAM_TO_HITS_local } from '../../../modules/local/bam_to_allsites_local/main_eshrem_single_end'
include { INSPIIRED_HITS_TO_SITESFINAL_local } from '../../../modules/local/insp_hits_to_sitesfinal_local/main'
include { ANNOTATE_SITES_local } from '../../../modules/local/annotate_sites_local/main_inspiired_until_alignment'

workflow POSTPROCESSING_twice_wfl {
    take:
    ch_aligned
    ch_processing_params
    ch_refGenome

    main:

    INDEX_SORT_BAM_local(ch_aligned)
    ch_sorted = INDEX_SORT_BAM_local.out.sorted

    //we need to merge the processing parameters with the aligned data
    ch_baminput = ch_sorted.combine(ch_processing_params, by: 0)

    BAM_TO_HITS_local(ch_baminput)
    ch_hits = BAM_TO_HITS_local.out.hits

    INSPIIRED_HITS_TO_SITESFINAL_local(ch_hits)
    ch_sitesfinal = INSPIIRED_HITS_TO_SITESFINAL_local.out.sitesfinal


    //we also need to join the channels to be able to give it to the process as input
    ch_annotinput = ch_sitesfinal.combine(ch_refGenome, by: 0)

    ch_sitesfinal.view { "Sitesfinal channel: ${it}" }
    ch_refGenome.view { "Reference genome channel: ${it}" }
    ch_annotinput.view { "Annotation input channel: ${it}" }

    ANNOTATE_SITES_local(ch_annotinput)
    ch_annotated = ANNOTATE_SITES_local.out.annotated

    emit:
    ch_annotated
}