# FINDVECTOR_local

## Description

When the CAR-T production has been made with Lentivirus, this process removes vector contamination from paired-end reads by identifying and filtering out read pairs where both R1 and R2 contain vector sequences.  
If the CAR-Ts have been produced with transposons, this process will be skipped.

---

## Script details

### 1. Loading libaries, functions and arguments

- The script uses the `Biostrings` and `ShortRead` libraries to perform the alignment and filtering.
- The script takes the following arguments:
  - `meta`: The sample name
  - `read1`: The path to the R1 FASTQ file
  - `read2`: The path to the R2 FASTQ file
  - `primerltr`: The primer and LTR bit of the vector concatenated
  - `vector`: The path to the vector FASTA file
  - `globalIdentity`: The minimum identity threshold for the alignment
- The script has two functions:
  - `calculateidentity`: This function calculates the identity of the alignment between the vector and the reads. Differently from INSPIIRED, this process uses minimap2 to perform the alignment due to its higher speed. However, minimap2 calculates the identity differently from BLAT, as it takes into account gaps to calculate the query length. Therefore, the identity is calculated as the number of matches divided by the query length.
  - `get_base_name`: This function extracts the base name from the header of the R1 and R2FASTQ files. The base name is the part before the first space, which read pairs have in common.

### 2. Aligning the reads

- The script uses minimap2 to align the reads to the vector.
- The script creates a data frame from the resulting PAF file from the alignment. The identity column of this dataframe is calculated as previously stated.
- The script filters the data frame to retain only the reads where the identity is greater than the minimum identity threshold.

### 3. Removing reads with vector in both reads

- As in INSPIIRED, the script merges the data frames for R1 and R2 to identify the reads that have vector in both reads.
- The script deletes the reads that have vector in both reads.

---

## Tools

| Tool        | Description                                                                 | Homepage                                                                 |
|------------|-----------------------------------------------------------------------------|-------------------------------------------------------------------------|
| minimap2   | Alignment tool for mapping DNA sequences                                    | [GitHub](https://github.com/lh3/minimap2)                               |
| biostrings | Biostrings is a set of low-level R functions for manipulating biological sequences | [Bioconductor](https://bioconductor.org/packages/release/bioc/html/Biostrings.html) |
| shortread  | ShortRead is a R package for input and manipulation of high-throughput sequence data | [Bioconductor](https://bioconductor.org/packages/release/bioc/html/ShortRead.html) |

---

## Input

| Name           | Type   | Description                           |
|----------------|--------|---------------------------------------|
| `sample`         | string | Sample name                            |
| `read1`          | file   | R1 FASTQ file                          |
| `read2`          | file   | R2 FASTQ file                          |
| `primer`         | string | Primer sequence                        |
| `ltrbit`         | string | LTR bit of the vector                  |
| `largeLTRFrag`   | string | Large LTR fragment flag                |
| `project_name`   | string | Project name                            |
| `mingDNA`        | string | Minimum genomic DNA threshold          |
| `vector_fasta`   | file   | Path to vector FASTA file               |

---

## Output

| Name   | Type | Description                     |
|--------|------|---------------------------------|
| `sample` | string | Sample name                     |
| `read1`  | file   | R1 FASTQ file with vector removed |
| `read2`  | file   | R2 FASTQ file with vector removed |

---

## Authors

- @liberentaizp
