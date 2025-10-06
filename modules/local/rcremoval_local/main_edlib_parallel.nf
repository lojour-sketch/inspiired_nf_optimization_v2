process RCremoval_edlib_parallel {

    publishDir '/home/lrenteria/inspiired_nf/results/9_reverse_complement_removal_regex', mode: 'symlink', overwrite: true

    cpus 8
    memory '60 GB'

    input:
    tuple val(meta), path(read1), path(read2), val(primer), val(ltrbit), val(largeLTRfrag), val(mingDNA)
    tuple val(sample), val(unique_linker), val(common_linker)

    output:
    tuple val(meta), path("*.rc_removed_R1.fastq.gz"), path("*.rc_removed_R2.fastq.gz"), emit: reads
    tuple val(meta), val(primer), val(ltrbit), val(largeLTRfrag), val(mingDNA)

    script:

    """
    #!/usr/bin/env python3

    import gzip
    import numpy as np
    from Bio.Seq import Seq
    import edlib
    from concurrent.futures import ThreadPoolExecutor
    import itertools
    import functools

        # Your parameters
    R1_in = "${read1}"
    R2_in = "${read2}"
    R1_out = "${meta}.rc_removed_R1.fastq.gz"
    R2_out = "${meta}.rc_removed_R2.fastq.gz"
    MAX_MISMATCH = 3
    MIN_LEN = ${mingDNA}

    

    def rc(seq):
        return str(Seq(seq).reverse_complement())
    
    largeLTRfrag_rc = rc("${largeLTRfrag}")
    common_linker_rc = rc("${common_linker}")

    def count_mismatches_cigar(cigar):
        #detects letters and numbers in cigar string and counts gap extensions
        total = 0
        num_str = ""
        
        for char in cigar:
            if char in '0123456789':
                num_str += char
            else:
                if num_str:
                    num = int(num_str)
                    # we only count gap extensions, as inspiired, so we count I, D and substract 1 from the length
                    if char in 'ID':
                        total += max(0, num - 1)
                    num_str = ""
        return total

    def process_read_batch(batch, marker_rc, MAX_MISMATCH, MIN_LEN):
        #Process a batch of reads
        results = []
        marker_len = len(marker_rc)
        half_marker = marker_len // 2
        
        for read_id, seq, qual in batch:
            # Fast edlib alignment
            aln = edlib.align(marker_rc, seq, mode="HW", task="path", k=MAX_MISMATCH)
            
            if not aln.get("locations"):
                # No alignment - trim from end
                cut = len(seq) - half_marker
                if cut >= MIN_LEN:
                    results.append((read_id, seq[:cut], qual[:cut]))
                continue
                
            start, end = aln["locations"][0]
            width = end - start + 1
            
            # Get CIGAR only when needed. if the alignment is none, we return what we did above
            if start > 0 or (start == 0 and width >= marker_len - 1):
                ## now we take the CIGAR part of the alignment to see how may mm, gaps... we have, and we can apply the INSPIIRED logic (matches+1, mm0, gaps0, gapextension-1)
                cigar_aln = edlib.align(marker_rc, seq, mode="HW", task="cigar", k=MAX_MISMATCH)
                mm = count_mismatches_cigar(cigar_aln["cigar"])
            else:
                mm = MAX_MISMATCH + 1  # we force a bad alignment because it does not meet inspiired criteria
                
            # we decide good and bad alignments based on INSPIIRED criteria
                # if the alignment starts after the first nt (start>0), good alignments need to have mm<=max_mismatch
                # if the alignment starts at the first nt (start=0), good alignments need to have mm<=max_mismatches and the aligned portion needs to be at least len(marker)-1 nt long
            good = (mm <= MAX_MISMATCH and start > 0) or (mm <= MAX_MISMATCH and start == 0 and width >= marker_len - 1)
            
            # we apply the cutting logic based on INSPIIRED criteria. 
                #if the alignment is good, we cut in the start of alignment
                #if the alignment is not good, we trim half of the length of the marker from the end of the read
            cut = start if good else len(seq) - half_marker
            
            #we make sure that the gDNA is longer than the mingDNA
            if cut >= MIN_LEN:
                results.append((read_id, seq[:cut], qual[:cut]))
        
        return results

    def fast_fastq_reader(filename):
        #Ultra-fast FASTQ reader
        with gzip.open(filename, 'rt') as f:
            while True:
                header = f.readline().strip()
                if not header: break
                seq = f.readline().strip()
                plus = f.readline().strip()
                qual = f.readline().strip()
                yield (header[1:], seq, qual)

    def main():     
        # Process files with thread pooling and batching
        def process_file(in_file, out_file, marker):
            with ThreadPoolExecutor(max_workers=8) as executor:
                batch_size = 10000
                
                with gzip.open(out_file, 'wt') as out_f:
                    reader = fast_fastq_reader(in_file)
                    
                    while True:
                        # Read batch without memory overload
                        batch = list(itertools.islice(reader, batch_size))
                        if not batch:
                            break
                        
                        # Process this batch
                        results = process_read_batch(batch, marker, MAX_MISMATCH, MIN_LEN)
                        for header, seq, qual in results:
                            out_f.write(f"@{header}\n{seq}\n+\n{qual}\n")
        
        # Process both files
        process_file(R1_in, R1_out, largeLTRfrag_rc)
        process_file(R2_in, R2_out, common_linker_rc)

    
    main()
    """
}