## LTRchecking_seqkit_local

### Description

This module checks if the R2 reads contain the primers and LTR sequences.  
We keep only the R2 reads that start with primer + LTR, and then we remove these sequences from all reads. In R1 reads, we remove all reads which pair did not contain the primer+LTR start.

**Script description:**

1. *Search for reads with primer + LTR*
   - We search for reads that START with primer+LTR in R2 reads with seqkit

2. *Remove paired reads from R1*
   - We remove all reads in R1 that pair with the R2 reads that did not have the primer+LTR start, to maintain the paired structure.
   - We do this by searching for read name

3. *Remove primer+LTR sequences from R2*
   - The R2 reads that we kept contain the primer+LTR sequence, and we need to remove it in otder to get the clean genomic DNA sequence.
   - We do this with seqkit subseq

---

### Tools

- seqkit
  - SeqKit is a cross-platform and ultrafast toolkit to process sequence files.
  - Homepage: https://bioinf.shenwei.me/seqkit/

---

### Input channel:

| Name        | Type    | Description |
|-------------|---------|-------------|
| `sample`      | string  | Sample name |
| `reads`       | list    | List that contains the paths to the R1 and R2 FASTQ read files. The first element is the R1 read file, and the second element is the R2 read file. |
| `primer`      | string  | Primer sequence |
| `LTR`         | string  | LTR sequence |
| `largeLTRFrag`| string  | Large LTR fragment sequence |
| `project`     | string  | Project name |
| `mingDNA`     | integer | Minimum genomic DNA length |

---

### Output channel: reads

| Name         | Type    | Description |
|--------------|---------|-------------|
| `sample`       | string  | Sample name |
| `R1`           | file    | Path to the cleaned R1 FASTQ read file |
| `R2`           | file    | Path to the cleaned R2 FASTQ read file |
| `primer`       | string  | Primer sequence |
| `LTR`          | string  | LTR sequence |
| `largeLTRFrag` | string  | Large LTR fragment sequence |
| `project`      | string  | Project name |
| `mingDNA`      | integer | Minimum genomic DNA length |

---

### Authors

- @liberentaizp
