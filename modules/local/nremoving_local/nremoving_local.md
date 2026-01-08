## N_REMOVING

### Description

This module removes the reads that contain N nucleotides.  
The posterior processes require reads without ambiguous N nucleotides, and as usually these type of reads are not too much, we will remove the reads that contain Ns.

**Script description:**

1. *Obtain read names of reads with Ns*
   - We obtain the read names of the reads that contain nts that are not A, C, G or T.
   - We do this for both R1 and R2 reads.

2. *Combine R1 and R2 headers*
   - As we want to remove the same reads from R1 and R2 to maintain the paired structure, we create a file that contains the read names of all R1 and R2 files with N nts.

3. *Remove reads with Ns and their pairs*
   - We convert the read names to R1 and R2 nomenclature
   - Then we remove the reads with Ns and their pairs from R1 and R2 files

---

### Tools

- seqkit
  - SeqKit is a cross-platform and ultrafast toolkit to process sequence files.
  - Homepage: https://bioinf.shenwei.me/seqkit/

---

### Input channel

| Name   | Type | Description |
|-------|------|-------------|
| `sample` | string | Sample name |
| `reads`  | list | List that contains the paths to the R1 and R2 FASTQ read files.

---

### Output channel: reads

| Name              | Type   | Description |
|-------------------|--------|-------------|
| `sample`            | string | Sample name |
| `N removed reads`   | file   | All the reads that have been filtered, without Ns.

---

### Authors

- @liberentaizp
