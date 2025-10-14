process RCremoval_local {

    publishDir '/home/lrenteria/inspiired_nf/results/9_reverse_complement_removal', mode: 'symlink', overwrite: true

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
    from Bio import SeqIO
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    import edlib
    from concurrent.futures import ProcessPoolExecutor
    from itertools import islice

    ######################### INPUTS #########################
    R1_in = "${read1}"
    R2_in = "${read2}"
    R1_out = "${meta}.rc_removed_R1.fastq.gz"
    R2_out = "${meta}.rc_removed_R2.fastq.gz"

    MAX_MISMATCH = 3
    MIN_LEN = ${mingDNA}   # e.g. 30 nt

    largeLTRfrag = "${largeLTRfrag}"
    common_linker = "${common_linker}"

    ######################### FUNCTIONS #########################

    # Reverse complement
    def rc(seq):
        return str(Seq(seq).reverse_complement())

    largeLTRfrag_rc = rc(largeLTRfrag)
    common_linker_rc = rc(common_linker)

    # Alignment + trimming logic
    def aligntrim(seq, qual, marker):
        # first we compute the location and cigar of the alignment with edlib to see where in the read it is. we will later need it
        aln = edlib.align(marker, seq, mode="HW", task="cigar", k=MAX_MISMATCH)
        
        if not aln["locations"]:  # no alignment at all means no good alignment
            cut = len(seq) - len(marker)//2
            return seq, qual

        start, end = aln["locations"][0]

        if start is None or end is None:
            cut = len(seq) - len(marker)//2
            return seq, qual


        width = end - start + 1

        # now we take the CIGAR part of the alignment to see how may mm, gaps... we have, and we can apply the INSPIIRED logic (matches+1, mm0, gaps0, gapextension-1)
        cigar = aln_cigar["cigar"]  # e.g., "10= 1I 20="

        #we count gaps only, and ignore mismatches, as in inspiired. 
        #we also ignore gap openings, that's why we add a -1 to each gap length
        mm = sum(max(0, int(op[:-1]) - 1) for op in cigar.split() if op[-1] in ("I", "D"))

        # we decide good and bad alignments based on INSPIIRED criteria
            # if the alignment starts after the first nt (start>0), good alignments need to have mm<=max_mismatch
            # if the alignment starts at the first nt (start=0), good alignments need to have mm<=max_mismatches and the aligned portion needs to be at least len(marker)-1 nt long
        good = (mm <= MAX_MISMATCH and start > 0) or (mm <= MAX_MISMATCH and start == 0 and width >= len(marker)-1)

        # we apply the cutting logic based on INSPIIRED criteria. 
            #if the alignment is good, we cut in the start of alignment
            #if the alignment is not good, we trim half of the length of the marker from the end of the read
        cut = start if good else len(seq) - len(marker)//2

        return seq[:cut], qual[:cut]


    # Process a chunk of SeqRecords
    def process_chunk(records, marker):
        trimmed = []
        for r in records:
            s, q = aligntrim(str(r.seq), r.letter_annotations["phred_quality"], marker)
            if len(s) >= MIN_LEN:
                new_r = SeqRecord(Seq(s), id=r.id, description=r.description,
                              letter_annotations={"phred_quality": q})
                trimmed.append(new_r)
        return trimmed

    # Yield FASTQ chunks
    def read_fastq_chunks(handle, chunk_size=100000):
        while True:
            chunk = list(islice(handle, chunk_size))
            if not chunk:
                break
            yield chunk

    #trimming function
    def trim_fastq_parallel(input_path, output_path, marker, cpus=8):
        with gzip.open(input_path, "rt") as handle, gzip.open(output_path, "wt") as out_handle:
            record_iterator = SeqIO.parse(handle, "fastq")
            # optional: parallel processing by chunks
            with ProcessPoolExecutor(max_workers=cpus) as executor:
                futures = {executor.submit(process_chunk, list(chunk), marker): i 
                        for i, chunk in enumerate(read_fastq_chunks(record_iterator))}
                
                # write chunks as they complete
                for fut in as_completed(futures):
                    trimmed_chunk = fut.result()
                    if trimmed_chunk:
                        SeqIO.write(trimmed_chunk, out_handle, "fastq")
    ######################### EXECUTION #########################

    # Run for paired-end
    trim_fastq_parallel(R1_in, R1_out, largeLTRfrag_rc, cpus=8)
    trim_fastq_parallel(R2_in, R2_out, common_linker_rc, cpus=8)

    """
}