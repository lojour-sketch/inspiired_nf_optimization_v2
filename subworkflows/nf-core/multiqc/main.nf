include { MULTIQC   } from '../../../modules/nf-core/multiqc/main'

workflow MULTIQC_wfl {
    take:
    fastqc_html
    fastqc_zip
    fastqc_html_trimmed
    fastqc_zip_trimmed
    multiqc_config
    extra_multiqc_config
    multiqc_logo
    replace_names
    sample_names

    main:
    // We want to join all the fastqc raw/trim html/zip files. they all have this format: [[id:xxxx, single_end=false], [paths/to/files]]
    
    //first we will join the raws in a channel
    ch_fastqc_raw = fastqc_html.join(fastqc_zip)
    
    //now we join the trimmed fastqcs in another channel
    ch_fastqc_trimmed = fastqc_html_trimmed.join(fastqc_zip_trimmed)
    
    //Now joining both
    ch_multiqc_flat = ch_fastqc_raw
        .join(
            ch_fastqc_trimmed
        )
        .map { meta, rawhtml, rawzip, trimhtml, trimzip -> [rawhtml, rawzip, trimhtml, trimzip].flatten() } // Strip metadata
        .flatten()                      //flatten nested lists to single paths



    //in order togive the process a channel with all the files at once
    ch_multiqc_collected = ch_multiqc_flat.collect()

    MULTIQC(
        ch_multiqc_collected,             //multiqc_files 
        multiqc_config,             //multiqc_config
        extra_multiqc_config,       //extra_multiqc_config
        multiqc_logo,               //multiqc_logo
        replace_names,              //replace_names
        sample_names                //sample_names
    )
    
    emit:
    report  = MULTIQC.out.report
    data    = MULTIQC.out.data
    plots   = MULTIQC.out.plots
    versions= MULTIQC.out.versions
}