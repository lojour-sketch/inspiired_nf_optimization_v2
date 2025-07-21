process LTRPRESENCE {
    input:
    # tuple val(meta), path(fastq), path(fastq_idx), path(undetermined), path(undetermined_idx), path(reports), path(stats), path(interop), path(versions)

    output:
    # tuple val(meta), path("LTRpresence_results/**.tsv"), emit: results

    script:
    """
    #!/usr/bin/env python3

    # First load samplesheet to get largeLTRFragment and primer sequences
    import pandas as pd

    samplesheet = pd.read_csv('${params.samplesheet}')
    




    """
}