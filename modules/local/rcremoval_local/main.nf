process RCremoval_local {

    publishDir '/../../../results/8_reverse_complement_removal', mode: 'copy', overwrite: true

    input:
    tuple val(sample), path(read1), path(read2)

    output:
    tuple val(sample), path("*.rc_removed_R1.fastq.gz"), path("*.rc_removed_R2.fastq.gz")

    script:
    def LTRbit = "GAAAATCTCTAGCA"
    def commonLinker = "CTCCGCTTAAGGGACT"
    //inside the awk process for each read we read the removing part from end to start and change each nt for their complement, then remove them from 
    """
    

    zcat ${read1} | awk -v seq="${LTRbit}" '
    BEGIN {
        rc=""
        for (i=length(seq); i>0; i--) {
            c=substr(seq,i,1)
            if(c=="A") c="T"; else if(c=="T") c="A"
            else if(c=="C") c="G"; else if(c=="G") c="C"
            else if(c=="a") c="t"; else if(c=="t") c="a"
            else if(c=="c") c="g"; else if(c=="g") c="c"
            rc=rc c
            }
        }
        {
            header = \$0; getline seq; getline plus; getline qual;
            gsub(rc, "", seq);
            qual = substr(qual, 1, length(seq));
            print header; print seq; print plus; print qual;
        }' | gzip > ${sample}.rc_removed_R1.fastq.gz
    
    zcat ${read2} | awk -v seq="${commonLinker}" '
    BEGIN {
        rc=""
        for (i=length(seq); i>0; i--) {
            c=substr(seq,i,1)
            if(c=="A") c="T"; else if(c=="T") c="A"
            else if (c=="C") c="G"; else if(c=="G") c="C"
            else if(c=="a") c="t"; else if(c=="t") c="a"
            else if(c=="c") c="g"; else if(c=="g") c="c"
            rc=rc c
            }
        }
       {
            header = \$0; getline seq; getline plus; getline qual;
            gsub(rc, "", seq);
            qual = substr(qual, 1, length(seq));
            print header; print seq; print plus; print qual;
        }' | gzip > ${sample}.rc_removed_R2.fastq.gz        

    """

}

//#rcLTRbit="\$(echo ${LTRbit} | tr ACGTacgt TGCAtgca | rev)"
//    #rccommonLinker="\$(echo ${rcLTRbit} | tr ACGTacgt TGCAtgca | rev)"

 //   #zcat ${read1} | awk 'NR%4==2 {gsub(/${rcLTRbit}/,"")} {print}' | gzip ${sample}.rc_removed_R1.fastq.gz

 //   #zcat ${read2} | awk 'NR%4==2 {gsub(/${rccommonLinker}/,"")} {print}' | gzip ${sample}.rc_removed_R2.fastq.gz
