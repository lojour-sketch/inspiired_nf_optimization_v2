#!/usr/bin/env python

#import modules
import sys
import pysam
import csv

#load the imputs
samfile = sys.argv[1]
allsites = sys.argv[2]
maxAlignStart = sys.argv[3]
minPctIdent = sys.argv[4]
maxFragLength = sys.argv[5]
sample = sys.argv[6]



#######################################################################################################################################################
#################################################################### QUALITY CHECK ####################################################################
#######################################################################################################################################################

# open the files
samFile = pysam.AlignmentFile(samfile, "rb")

num_reads = samFile.count(until_eof=True)
print("Number of reads in BAM:", num_reads)

# Close and reopen to reset iterator
samFile.close()

# reopen the files
samFile = pysam.AlignmentFile(samfile, "rb")


with open(f"{sample}_tmpFile.tsv", "w", newline="") as outfile:
    writer = csv.writer(outfile, delimiter="\t")
    # write header
    writer.writerow(["readname","chr","from","strand","start", "end", "cigar","qStart","PercIdent","tBaseInsert","flag", "tags"])
    print("Entering for loop in: ", f"{sample}_tmpFile.tsv")
    # check if all reads have >minpercIdent, <qStart and <tBasesInsert, and are properly paired
    for read in samFile.fetch(until_eof=True):
        if read.is_read1:
            read_num = "R1"
        elif read.is_read2:
            read_num = "R2"
        if read.is_reverse:
            strand = "-"
        elif not read.is_reverse:
            strand = "+"
        #cigar tuples save the cigar string as a list of tuples, where the first element is the operation and the second is the length
        # first we will calculate the qStart, which is the equivalent pf the initial soft-clippiong bases in cigar
        cigartuple = read.cigartuples
        qStart = 0
        totalMatches = 0
        tBaseInsert = 0
        if cigartuple[0][0] == 4:
            qStart = cigartuple[0][1]
        for cigar in cigartuple:
            if cigar[0] == 0:
                totalMatches += cigar[1]
            elif cigar[0] == 2:
                #tBaseInsert variable in blat alignments means the number of bases inserted in the target (reference)
                #however, in CIGAR nomenclature, these are represented as deletions in the query (read), so we take the cigartuples deletion operation as tBaseInsert
                tBaseInsert += cigar[1]
        queryLength = read.query_length
        if queryLength==0 or queryLength==None: #because some reads may have 0 length, we avoid division by zero. also, pysam may return None for unmapped reads
            queryLength = 1
        percIdent = 100 * totalMatches / queryLength
        # check if the read is mapped , paired, not duplicated, not secondary, not supplementary, has a minimum percent identity and inserts, and has a maximum fragment length
            #supplementary reads: reads that are aligned by parts, this is, reads that have been split into multiple parts
            #secondary reads: reads that align to more than one site
        if read.is_proper_pair and read.is_paired and not read.is_duplicate and not read.is_secondary and not read.is_supplementary and float(percIdent) > float(minPctIdent) and int(tBaseInsert) < 5 and read.template_length < int(maxFragLength) and int(qStart) < int(maxAlignStart):
            writer.writerow([read.query_name, read.reference_name, read_num, strand, read.reference_start, read.reference_end, read.cigarstring, qStart, percIdent, tBaseInsert, read.flag, read.get_tags()])

#######################################################################################################################################################
############################################################### POTENTIAL INSERTION SITES #############################################################
#######################################################################################################################################################

#we already have a tsv files with all the filtered reads. now we have to merge the reads that correspond to the same pair. 
#as the bam file is sorted by read name, we can iterate through the file and merge the reads with the same name

with open(f"{sample}_tmpFile.tsv", "r", newline="") as infile, open(f"{sample}_allSites.tsv", "w", newline="") as outfile:
    reader = csv.DictReader(infile, delimiter="\t")
    writer = csv.writer(outfile, delimiter="\t")
    #write header
    writer.writerow(["seqnames", "start", "end", "strand", "revmap", "pairingID","samplename", "ID"])

    last_readname = None
    pair_reads = []
    pairingID = 0

    for i, row in enumerate(reader, start=1):
        #we save the line number to later use it in revmap
        row["line_num"] = i
        readname = row["readname"]

        if last_readname is None:
            last_readname = readname
        
        if readname != last_readname:
            if len(pair_reads) == 2:
                #we have a pair of reads, so we can merge them.
                #we have to assign the reads to r1 and r2 based on the read_num, 
                #make sure that their line_num is different (different reads) 
                #and that their readname is the same (pairs)
                r1 = [r for r in pair_reads if r["from"] == "R1"][0]
                r2 = [r for r in pair_reads if r["from"] == "R2"][0]
                if r1["line_num"] != r2["line_num"] and r1["readname"] == r2["readname"]:
                    chrname = r1["chr"]
                    strand = r2["strand"]

                    start_r1 = int(r1["start"])
                    start_r2 = int(r2["start"])
                    #We will assign the start position to the R2 start when the strand is +, and to the R1 start when the strand is -
                    #we do this because we cannot create a GRanges object with start > end, so we will always have start < end, but we will keep the strand information
                    start = start_r2 if strand == "+" else start_r1
                    end = start_r1 if strand == "+" else start_r2
                    
                    #to create the revmap, we use the line numbers of the reads in the tmp file
                    revmap = f"{pair_reads[0]['line_num']},{pair_reads[1]['line_num']}"
                    ranges = f"{start}-{end}"
                    pairingID += 1
                    writer.writerow([chrname, start, end, strand, revmap, pairingID, sample, last_readname])

            #reset for next read
            pair_reads = []
            last_readname = readname
        #this is for the last read. we need a 'next read' to trigger the write, so we do it after the loop
        pair_reads.append(row)
    if len(pair_reads) == 2:
        r1 = [r for r in pair_reads if r["from"] =="R1"][0]
        r2 = [r for r in pair_reads if r["from"] =="R2"][0]
        chrname = r1["chr"]
        strand = r2["strand"]

        start_r1 = int(r1["start"])
        start_r2 = int(r2["start"])
        #We will assign the start position to the R2 start when the strand is +, and to the R1 start when the strand is -
        #we do this because we cannot create a GRanges object with start > end, so we will always have start < end, but we will keep the strand information
        start = start_r2 if strand == "+" else start_r1
        end = start_r1 if strand == "+" else start_r2

        revmap = f"{pair_reads[0]['line_num']},{pair_reads[1]['line_num']}"
        ranges = f"{start}-{end}"
        pairingID += 1

        writer.writerow([chrname, start, end, strand, revmap, pairingID, sample, last_readname])

#close the sam file
samFile.close()
