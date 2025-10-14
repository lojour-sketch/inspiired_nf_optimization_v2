
// first draft?

// workflow {
//     main:
    
//     // Process 1: Quality Filter SAM files
//     QUALITY_FILTER_SAM(params.r1_sam, params.r2_sam)
    
//     // Process 2: Pair filtered alignments  
//     READ_PAIRING(QUALITY_FILTER_SAM.out.filtered_sam, params.keys_file)
    
//     // Process 3: Classify pairs
//     CLASSIFY_PAIRS(READ_PAIRING.out.valid_pairs)
    
//     // Process 4: Generate outputs
//     GENERATE_OUTPUTS(CLASSIFY_PAIRS.out.classified)
// }