#!/usr/bin/env python

#import modules
import sys
import pysam
import csv

#load the imputs
samfile_r1 = sys.argv[1]
samfile_r2 = sys.argv[2]
maxAlignStart = sys.argv[3]
minPctIdent = sys.argv[4]
maxFragLength = sys.argv[5]
sample = sys.argv[6]



#######################################################################################################################################################
#################################################################### QUALITY CHECK ####################################################################
#######################################################################################################################################################

# open the files
samFile = pysam.AlignmentFile(samfile_r1, "rb")

num_reads = samFile.count(until_eof=True)
print("Number of reads in R1 BAM:", num_reads)

# Close and reopen to reset iterator
samFile.close()

# reopen the files
samFile = pysam.AlignmentFile(samfile_r1, "rb")


with open(f"{sample}_R1_hits.tsv", "w", newline="") as outfile:
    writer = csv.writer(outfile, delimiter="\t")
    # write header
    writer.writerow(["read_id","seqnames","start","end","width", "strand", "matches","qStart", "qEnd", "qSize", "qName", "tBaseInsert", "from"])
    print("Entering for loop in: ", f"{sample}_tmpFile.tsv")
    
    for read in samFile.fetch(until_eof=True):
        read_num = "R1"
        #cigar tuples save the cigar string as a list of tuples, where the first element is the operation and the second is the length
        # first we will calculate the qStart, which is the equivalent pf the initial soft-clippiong bases in cigar
        cigartuple = read.cigartuples
        qStart = 0
        totalMatches = 0
        tBaseInsert = 0
        queryLength = read.query_length
        if cigartuple[0][0] == 4:
            qStart = cigartuple[0][1]
        else:
            qStart = 0
        if cigartuple[-1][0] == 4:
            qEnd = cigartuple[-1][1]
        else:
            qEnd = queryLength
        for cigar in cigartuple:
            if cigar[0] == 0:
                totalMatches += cigar[1]
            elif cigar[0] == 2:
                #tBaseInsert variable in blat alignments means the number of bases inserted in the target (reference)
                #however, in CIGAR nomenclature, these are represented as deletions in the query (read), so we take the cigartuples deletion operation as tBaseInsert
                tBaseInsert += cigar[1]
        #now we save the strand and also define start and end depending on the processBLATData in INSPIIRED
        tStart = read.reference_start
        tEnd = read.reference_end
        if read.is_reverse:
            strand = "-"
            start = tStart - (queryLength - qEnd - 1)
            end = tEnd + qStart
        elif not read.is_reverse:
            strand = "+"
            start = tStart - qStart
            end = tEnd + (queryLength - qEnd - 1)
        if queryLength==0 or queryLength==None: #because some reads may have 0 length, we avoid division by zero. also, pysam may return None for unmapped reads
            queryLength = 1
        percIdent = 100 * totalMatches / queryLength
        width = end - start
        if not read.is_secondary and not read.is_supplementary and float(percIdent) > float(minPctIdent) and int(tBaseInsert) < 5 and read.template_length < int(maxFragLength) and int(qStart) < int(maxAlignStart):
            writer.writerow([read.query_name, read.reference_name, start, end, width, strand, totalMatches, qStart, qEnd,queryLength, read.query_name, tBaseInsert, read_num])

samFile.close()

# now for R2

# open the files
samFile = pysam.AlignmentFile(samfile_r2, "rb")

num_reads = samFile.count(until_eof=True)
print("Number of reads in R2 BAM:", num_reads)

# Close and reopen to reset iterator
samFile.close()

# reopen the files
samFile = pysam.AlignmentFile(samfile_r2, "rb")


with open(f"{sample}_R2_hits.tsv", "w", newline="") as outfile:
    writer = csv.writer(outfile, delimiter="\t")
    # write header
    writer.writerow(["read_id","seqnames","start","end","width", "strand", "matches","qStart", "qEnd", "qSize", "qName", "tBaseInsert", "from"])
    print("Entering for loop in: ", f"{sample}_tmpFile.tsv")
    
    for read in samFile.fetch(until_eof=True):
        read_num = "R2"
        #cigar tuples save the cigar string as a list of tuples, where the first element is the operation and the second is the length
        # first we will calculate the qStart, which is the equivalent pf the initial soft-clippiong bases in cigar
        cigartuple = read.cigartuples
        qStart = 0
        totalMatches = 0
        tBaseInsert = 0
        queryLength = read.query_length
        if cigartuple[0][0] == 4:
            qStart = cigartuple[0][1]
        else:
            qStart = 0
        if cigartuple[-1][0] == 4:
            qEnd = cigartuple[-1][1]
        else:
            qEnd = queryLength
        for cigar in cigartuple:
            if cigar[0] == 0:
                totalMatches += cigar[1]
            elif cigar[0] == 2:
                #tBaseInsert variable in blat alignments means the number of bases inserted in the target (reference)
                #however, in CIGAR nomenclature, these are represented as deletions in the query (read), so we take the cigartuples deletion operation as tBaseInsert
                tBaseInsert += cigar[1]
        #now we save the strand and also define start and end depending on the processBLATData in INSPIIRED
        tStart = read.reference_start
        tEnd = read.reference_end
        if read.is_reverse:
            strand = "-"
            start = tStart - (queryLength - qEnd - 1)
            end = tEnd + qStart
        elif not read.is_reverse:
            strand = "+"
            start = tStart - qStart
            end = tEnd + (queryLength - qEnd - 1)
        if queryLength==0 or queryLength==None: #because some reads may have 0 length, we avoid division by zero. also, pysam may return None for unmapped reads
            queryLength = 1
        percIdent = 100 * totalMatches / queryLength
        width = end - start
        if not read.is_secondary and not read.is_supplementary and float(percIdent) > float(minPctIdent) and int(tBaseInsert) < 5 and read.template_length < int(maxFragLength) and int(qStart) < int(maxAlignStart):
            writer.writerow([read.query_name, read.reference_name, start, end, width, strand, totalMatches, qStart, qEnd,queryLength, read.query_name, tBaseInsert, read_num])

samFile.close()
