 process N_REMOVING {

    tag "${meta.id}"
    publishDir "${params.runfolderDir}/../results/5_removed_n", mode: 'copy', overwrite: true

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*_N.fq.gz"), emit: reads

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """

        # Extract headers of reads containing symbols different from ACGT in both files
        echo "Number of reads with strange symbols in R1:"
        zcat ${prefix}_1_val_1.fq.gz | awk 'NR%4==1 {sub(/^@/,""); header=\$0} NR%4==2 {seq=\$0} NR%4==0 {if(seq !~ /^[ACGT]*\$/) print header}' > R1_N_headers.txt
        

        echo "Number of reads with strange symbols in R2:"
        zcat ${prefix}_2_val_2.fq.gz | awk 'NR%4==1 {sub(/^@/,""); header=\$0} NR%4==2 {seq=\$0} NR%4==0 {if(seq !~ /^[ACGT]*\$/) print header}' > R2_N_headers.txt
        


        # Combine R1 and R2 headers for each output
        cat R1_N_headers.txt R2_N_headers.txt | sort | uniq > N_headers_all.txt


        # to find the headers in R1 we will change the number in {num}:N:0 to 1
        sed 's/[12]:N:0/1:N:0/' N_headers_all.txt > N_headers_clean_R1.txt
        sed 's/[12]:N:0/2:N:0/' N_headers_all.txt > N_headers_clean_R2.txt


        # Remove all reads with N from both R1 and R2
        seqkit grep -n -v -f N_headers_clean_R1.txt ${prefix}_1_val_1.fq.gz | gzip > ${prefix}_1_N.fq.gz
        seqkit grep -n -v -f N_headers_clean_R2.txt ${prefix}_2_val_2.fq.gz | gzip > ${prefix}_2_N.fq.gz


    """
 }
    
    
    
    
    
    
    
    
    