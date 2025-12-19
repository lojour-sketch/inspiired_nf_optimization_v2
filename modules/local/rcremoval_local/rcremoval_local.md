## RCremoval_inspiired_local

### Description

This process removes the reverse complement of the *other read's* common linker or largeLTRFragment.

**Justification of the process:**

When R1 and R2 reads are sequenced, it may happen that (if the sequencing fragment is short) the read1 sequencing arrives until the read2 large LTR fragment, sequencing the virus genome along with the read1 sequence. Likewise, if the sequencig fragment is short, the read2 sequencing can arrive until the read1 primers.  
As we only want to keep the clean genomic DNA, we have to remove these artificial fragments. For that, we will search for the reverse complement of the large LTR fragment in the read1 reads and remove them (with all that comes after). Correspondingly we will also search for the reverse complement of the read1 common linker in read2 reads, and remove them (with all that comes after).

**Process' script description:**

1. *See if files are empty*
   - As sometimes control files (that are frequently empty) are present in the pipeline, we will first check if the files are empty.
   - If they are empty, we will skip the process and just create empty output files.
   - If they are not empty, we will continue with the process, running the R script with the following workflow

**R script description:**

1. *Load libraries and arguments*
   - Libraries are: Biostrings, ShortRead, pwalign
   - Arguments are: meta, read1, read2, largeLTRfrag, common_linker

2. *Check if files are empty*
   - We already check this in the process' script, but we will do it again here to be sure
   - If the files are empty, we will create empty output files

3. *Verify that reads are paried*
   - We verify that the R1 and R2 reads have the same length, and if not, we stop the process

4. *Alignment and obtaining cut point for each read*
   - We obtain the reverse complements of the markers
   - We first make a function that converts the output of the pairwise alignment into a dataframe
   - We then apply our function to perform the alignment and apply the INSPIIRED criteria
     - We assume no reads have good alignments
     - Ensure start position is at least 1
     - We apply INSPIIRED criteria:
       - IF the overlap is somewhere NOT in the start of the read, we need the alignment to have mm < maxMismatch to be considered as good
       - Now if the overlap is in the start of the read, we need mm < masMismatch AND the aligned portion has to be longer than length(marker)-1 to be considered as good
     - We assign nchar(reads) - nchar(marker)/2 to all reads first, and then we change the value for the good alignments
       - We do this to precautionarily trim half the length of the marker from the final of the reads with no good alignment, to avoid having incomplete markers at the end that would not give a good alignment but need to be removed
     - For the alignments considered as good, we trim from the start of the alignment. This way we remove the aligned part and all that comes after
     - We ensure that the cut position is at least 1 and that it does not exceed read length

5. *Trim reads*
   - We trim the reads from the cut_vector point
   - For R1 reads we find the reverse complement of the largeLTRfrag
   - For R2 reads we find the reverse complement of the common_linker
   - We also trim the quality sequences the same way as the sequence

6. *Verify synchronization*
   - We verify that the length of both R1 and R2 reads is the same, and we create the output files.

---

### Tools

- Biostrings
  - Biostrings is a set of low-level R functions for manipulating biological strings.
  - Homepage: https://bioconductor.org/packages/release/bioc/html/Biostrings.html
- ShortRead
  - ShortRead is a package for input, output, manipulation, and quality control of high-throughput sequence data.
  - Homepage: https://bioconductor.org/packages/release/bioc/html/ShortRead.html
- pwalign
  - pwalign is a package for pairwise sequence alignment.
  - Homepage: https://bioconductor.org/packages/release/bioc/html/pwalign.html

---

### Input channel

| Name           | Type    | Description |
|----------------|---------|-------------|
| `sample`         | string  | Sample name |
| `read1`          | file    | Path to the R1 FASTQ read file |
| `read2`          | file    | Path to the R2 FASTQ read file |
| `primer`         | string  | Primer sequence |
| `ltrbit`         | string  | LTR bit sequence |
| `largeltrfrag`   | string  | Large LTR fragment sequence |
| `project`        | string  | Project name |
| `mingDNA`        | integer | Minimum genomic DNA length |
| `unique_linker`  | string  | Unique linker sequence |
| `common_linker`  | string  | Common linker sequence |

---

### Output channel: reads

| Name   | Type | Description |
|--------|------|-------------|
| `sample` | string | Sample name |
| `R1`     | file   | Path to the cleaned R1 FASTQ read file |
| `R2`     | file   | Path to the cleaned R2 FASTQ read file |

---

### Authors

- @liberentaizp
