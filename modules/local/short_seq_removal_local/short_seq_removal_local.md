## SHORTREMOVE_local

### Description

This module removes the reads that are shorter than the minimum length defined in the samplesheet.  
We also remove the pairs of the short reads, to maintain the paired structure.

**Script description:**

1. *Remove reads shorter than the minimum length*
   - We use seqkit seq -m to remove the too short sequences  
     - https://bioinf.shenwei.me/seqkit/usage/#seq

2. *Check the file stats: size and content*
   - We check the number of reads in the files with seqkit stats  
     - https://bioinf.shenwei.me/seqkit/usage/#stats
   - We check if the files are empty after the filtering. if they are, we create empty output files

3. *Remove pairs of the short reads*
   - To maintain the paired structure, we remove the pairs of the too short reads  
     - https://bioinf.shenwei.me/seqkit/usage/#pair
   - We check if the pairing was succesful. If not, we stop the process

---

### Tools

- **seqkit**  
  SeqKit is a cross-platform and ultrafast toolkit to process sequence files.  
  https://bioinf.shenwei.me/seqkit/

---

### Input channel

| Name         | Type    | Description |
|--------------|---------|-------------|
| `sample`       | string  | Sample name |
| `read1`        | file    | Path to the R1 FASTQ read file |
| `read2`        | file    | Path to the R2 FASTQ read file |
| `primer`       | string  | Primer sequence |
| `ltrbit`       | string  | LTR bit sequence |
| `largeltrfrag` | string  | Large LTR fragment sequence |
| `project`      | string  | Project name |
| `mingDNA`      | integer | Minimum genomic DNA length |

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
