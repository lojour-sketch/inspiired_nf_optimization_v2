# INSPIIRED modules

## Modules

## normalize_index_length_local

### Description

In various INSPIIRED runs, some sample unique linkers have different lengths. We need this sequence to demultiplex the data by samples, and we need it to be the same length for all samples.  
As the iSL (linker sequence) in INSPIIRED runs is like:  
`<Sample unique linker><UMI sequence><Common linker>`

For the samples that have a longer sample sequence, we will remove the extra nucleotides and add them to the UMI sequence. For this, we remove the nts from the longest sample linkers and we save the names of their corresponding samples to later expand the length of their UMIs.

**Script description:**

1. *Load libraries and arguments*
   - Libraries are: csv, sys, argparse, pathlib
   - Arguments are: samplesheet

2. *Process index lengths*
   - By parsing the samplesheet, we get the lengths of all the sample unique indexes (index column in samplesheet).
   - We save the names of the samples that contained longer indexes
   - We remove one nt from the longest indexes.

3. *Rewrite samplesheet*
   - We rewrite the samplesheet with the new index lengths

---

### Input channel

| Name        | Type   | Description |
|-------------|--------|-------------|
| `sample`      | string | Sample name |
| `samplesheet` | file   | Path to samplesheet CSV |

---

### Output channel: normalized

| Name                  | Type   | Description |
|-----------------------|--------|-------------|
| `sample`               | string | Sample name |
| `normalized_samplesheet`| file   | Path to the normalized samplesheet CSV |
| `runfolder`            | path   | Path to the runfolder |

### Output channel: modified_samples

| Name             | Type | Description |
|------------------|------|-------------|
| `modified_samples` | file | Optional txt file containing the names of the samples that had their index length modified.

---


# CREATE_demux_samplesheet_local

**Description:**  
Creates a samplesheet for demultiplexing when initial inputs are Undetermined FASTQ files.

---

## Script Details

The module runs different scripts depending on the sequencing instrument, following Illumina guidelines for reverse complementing i5 sequences:  

- **MiSeq:** index2 sequence must be forward-oriented  
- **NextSeq:** index2 sequence must be reverse complemented

<details>
<summary>1. Forward script</summary>

- Creates a samplesheet with columns: `sample_id` and `barcode`  
  - `barcode` = concatenation of index1 + index2  
- Index2 sequence is in the same orientation as index1

</details>

<details>
<summary>2. Reverse script</summary>

- Creates a samplesheet with columns: `sample_id` and `barcode`  
  - `barcode` = concatenation of index2 + reverse complement of index1  
- Index2 sequence is reverse complemented

</details>

---

## Tools

<details>
<summary>python</summary>

**Description:** Python interpreter  
**Homepage:** [https://www.python.org/](https://www.python.org/)

</details>

---

## Input

| Name | Type | Description | Pattern |
|------|------|-------------|---------|
| `sample` | string | Sample name | — |
| `normalized_samplesheet` | file | Normalized samplesheet with indexes of same length across all samples | — |
| `rundir` | directory | Run directory | — |

---

## Output

| Name | Type | Description | Pattern |
|------|------|-------------|---------|
| `demux_sheet` | file | Demultiplexing samplesheet used by `fastqk` to demultiplex reads | `DemuxSampleSheet.tsv` |

---

# CREATE_demux_samplesheet_local

**Description:**  
Creates a samplesheet for demultiplexing when initial inputs are Undetermined FASTQ files.

---

## Script Details

The module runs different scripts depending on the sequencing instrument, following Illumina guidelines for reverse complementing i5 sequences:  

- **MiSeq:** index2 sequence must be forward-oriented  
- **NextSeq:** index2 sequence must be reverse complemented

<details>
<summary>1. Forward script</summary>

- Creates a samplesheet with columns: `sample_id` and `barcode`  
  - `barcode` = concatenation of index1 + index2  
- Index2 sequence is in the same orientation as index1

</details>

<details>
<summary>2. Reverse script</summary>

- Creates a samplesheet with columns: `sample_id` and `barcode`  
  - `barcode` = concatenation of index2 + reverse complement of index1  
- Index2 sequence is reverse complemented

</details>

---

## Tools

<details>
<summary>python</summary>

**Description:** Python interpreter  
**Homepage:** [https://www.python.org/](https://www.python.org/)

</details>

---

## Input

| Name | Type | Description | Pattern |
|------|------|-------------|---------|
| `sample` | string | Sample name | — |
| `normalized_samplesheet` | file | Normalized samplesheet with indexes of same length across all samples | — |
| `rundir` | directory | Run directory | — |

---

## Output

| Name | Type | Description | Pattern |
|------|------|-------------|---------|
| `demux_sheet` | file | Demultiplexing samplesheet used by `fastqk` to demultiplex reads | `DemuxSampleSheet.tsv` |

---

# BCL2FASTQ_local

**Description:**  
Converts Illumina BCL files to FASTQ format with demultiplexing.  
Performs base calling and demultiplexing of Illumina sequencing runs using `bcl2fastq2`, generating sample-specific FASTQ files based on index sequences defined in the sample sheet.

---

## Parameter Configuration and Processing

1. Use-bases-mask

- Custom bases mask for read structure: `I20Y159,I12,Y143`  
  - Read 1: 20bp index (I20) + 159bp data (Y159)  
  - Index 1: 12bp index (I12)  
  - Read 2: 143bp data (Y143)

2. Barcode mismatches

- Allowed mismatches: `2,2`  
  - 2 mismatches in first index  
  - 2 mismatches in second index  
  - Ensures higher read recovery during demultiplexing

3. No lane splitting

- Single line used for sequencing (no lane splitting)

4. File organization
- Creates sample-specific subdirectories:  
  `results/project/sample/fastq.gz`

5. Quality reports

- Reports/: HTML per-sample metrics  
- Stats/: JSON detailed statistics  
- InterOp/: Binary files for Illumina SAV viewer

**Important notes:**
- Independently of the sequencing instrument, BCL accepts always forward oriented barcodes. Different instruments require forward or reverse complement i5 indexes (Golay sequences in our case), but BCL2FASTQ handles this automatically.

<details>
<summary>Performance optimization</summary>

- `-r 25`: threads for reading BCL  
- `-p 25`: threads for processing  
- `-w 25`: threads for writing FASTQ  
- Total CPU usage ~75 cores at peak

</details>

---

## Tools

<details>
<summary>bcl2fastq</summary>

**Description:** Illumina BCL to FASTQ converter  
**Homepage:** [https://support.illumina.com/sequencing/sequencing_software/bcl2fastq-conversion-software.html](https://support.illumina.com/sequencing/sequencing_software/bcl2fastq-conversion-software.html)

</details>

---

## Input

| Name | Type | Description | Pattern |
|------|------|-------------|---------|
| `sample` | string | Sample name | — |
| `primer` | string | Primer information from SampleSheet | — |
| `ltrbit` | string | LTR bit info from SampleSheet | — |
| `largeLTRFrag` | string | Large LTR fragment flag from SampleSheet | — |
| `project_name` | string | Project name for output organization | — |
| `mingDNA` | string | Minimum genomic DNA threshold | — |
| `meta` | map | Metadata map containing additional run information | — |
| `samplesheet` | file | Illumina CSV sample sheet with required columns: Sample_ID, index, index2, common_linker, primer, ltrbit, largeLTRFrag, Sample_Project, mingDNA, minPctIdent, maxAlignStart, maxFragLength, refGenome, vectorSeq | `*.csv` |
| `run_folder` | directory | Illumina run folder containing BCL files in required structure (`Data/Intensities/BaseCalls/L00x`, `InterOp/`, `RunInfo.xml`, `RunParameters.xml`) | — |

---

## Output

| Name | Type | Description | Pattern |
|------|------|-------------|---------|
| `meta` | map | Metadata map propagated from input | — |
| `sample_fastq` | file | Demultiplexed FASTQ per sample. `{Sample_Name}_S{Sample_Number}_R{Read}_001.fastq.gz` | `*_R{1,2}_001.fastq.gz` |
| `index_fastq` | file | Index read FASTQ files (optional). `{Sample_Name}_S{Sample_Number}_I{Index}_001.fastq.gz` | `*_I{1,2}_001.fastq.gz` |
| `undetermined_fastq` | file | Reads that could not be assigned to samples. `Undetermined_S0_R{1,2}_001.fastq.gz` | `Undetermined_S0_R{1,2}_001.fastq.gz` |
| `undetermined_index` | file | Index sequences for undetermined reads. `Undetermined_S0_I{1,2}_001.fastq.gz` | `Undetermined_S0_I{1,2}_001.fastq.gz` |
| `reports_dir` | directory | HTML demultiplexing reports | `Reports/` |
| `stats_dir` | directory | JSON detailed statistics | `Stats/` |
| `interop_files` | file | Binary InterOp files for Illumina SAV viewer | `InterOp/*.bin` |

---

## UMI_EXTRACT_local

### Description

This module extracts the UMI sequence from the read1 FASTQ file.  
As we removed the sample linker from all read1s in the previous process, the UMI sequence is the first 12 (or 13 if we removed one nt from the sample linker) of the R1 reads. We save the UMI sequence in the read headers.  

We need to remove all reads that do not contain the `<UMI><common_linker>` sequence in the first part. Using `umi_tools extract` we:  
1. Remove all reads that do not contain the `<UMI><common_linker>` sequence in the first part  
2. For reads that contain both sequences, save the UMI sequence in the read header and remove both sequences from the read, leaving the clean genomic DNA start intact  

**Parameter configuration:**

1. *Extraction method:*  
   - Use `regex` option, as the barcode pattern contains regular expressions  

2. *Bc-pattern:*  
   - Pattern: `(?P<umi_1>[ATCGN]{${umi_length}})(?P<cell_2>${linker2})`  
   - Ensures all reads contain the UMI sequence and common linker; reads that do not match are removed  

3. *Read inputs and outputs:*  
   - `I`: Read1 input file  
   - `S`: Read1 output file  
   - `Read2-in`: Read2 input file  
   - `Read2-out`: Read2 output file  

4. *Log file:*  
   - Saved in the same directory as the output files

---

### Tools

| Tool | Description | Homepage |
|------|-------------|---------|
| umi_tools | Tool for handling Unique Molecular Identifiers (UMIs) from 10X Genomics data | [link](https://umi-tools.readthedocs.io/en/latest/reference/extract.html?highlight=extract#module-umi_tools.extract) |

---

### Input channel

| Name | Type | Description |
|------|------|-------------|
| `sample_id` | string | Sample name |
| `linker1` | string | Sample unique linker |
| `linker2` | string | Common linker |
| `reads` | list | List of paths to R1 and R2 FASTQ read files. First element: R1, second: R2 |
| `was_modified` | boolean | True if the sample unique linker was modified, False otherwise |

---

### Output channel

| Name | Type | Description |
|------|------|-------------|
| `sample_id` | string | Sample name |
| `R1` | file | Path to the UMI extracted R1 FASTQ read file |
| `R2` | file | Path to the UMI extracted R2 FASTQ read file |
| `log` | file | Path to the log file |

---

(FASTQC here, already documented in the internet)


---

(TRIMGALORE here, already documented in the internet)

---

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

# GENOME_INDEXING

## Description

This process indexes the reference genome using bowtie2, to be able to do the alignment a posteriori.

---

### Parameter configuration

- **runMode**: This is the run mode of STAR. This time we use the `--runMode genomeGenerate` option, which creates the indexed genome.  
- **genomeDir**: Directory that will contain the indexed genome.  
- **genomeFastaFiles**: Path to the reference genome fasta file.  
- **runThreadN**: Number of threads to use for indexing.  

---

## Tools

| Tool | Description                                   | Homepage |
|------|-----------------------------------------------|----------|
| STAR | Alignment tool for mapping DNA sequences     | [GitHub](https://github.com/alexdobin/STAR) |

---

## Input

| Name           | Type   | Description                       |
|----------------|--------|-----------------------------------|
| `genome_name`    | string | Name of the genome                |
| `refGemomeFile`  | file   | Path to the reference genome fasta file |

---

## Output

| Name              | Type      | Description                                      |
|------------------|-----------|--------------------------------------------------|
| `genome_name`       | string    | Name of the genome                               |
| `indexedGenomeDir`  | directory | Directory that will contain the indexed genome  |

---

# ALIGNMENT_local

**Description:**  
Aligns the R1 and R2 reads to the reference genome in a **paired-end manner**, to identify genomic regions where insertions occurred.  
INSPIIRED aligns both reads separately, then merges them after filtering. Aligning in pairs reduces multihits and chimeras.  
This module uses the STAR aligner for mapping.

---

## Tools

<details>
<summary>STAR</summary>

**Description:**  
STAR is a fast universal RNA-seq aligner. In this module, it aligns R1 and R2 reads to the reference genome to locate insertion sites.

**Homepage:** [https://github.com/alexdobin/STAR](https://github.com/alexdobin/STAR)  
**Documentation:** [STAR manual](https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf)  
**DOI:** doi.org/10.1093/bioinformatics/bts635

</details>

---

## Input

| Name | Type | Description |
|------|------|-------------|
| `sample` | string | Sample name |
| `r1` | file | FastQ file containing all R1 reads |
| `r2` | file | FastQ file containing all R2 reads |
| `genome_index` | string | Path to the genome index generated by `GENOME_INDEXING_local` process |

---

## Output

| Name | Type | Description | Pattern |
|------|------|-------------|---------|
| `sample` | string | Sample name | — |
| `*.Aligned.*` | file | STAR alignment output files (BAM, logs, etc) | `*.Aligned.*` |

*Tip:* Use backticks (\``) to highlight file patterns or code elements.

---

# INDEX_SORT_BAM_local

## Description

This process sorts the BAM file by coordinate and indexes it with samtools, in order to be able to perform the processing a posteriori.

---

### Command line structure

- First we index the BAM file with `samtools index`. For this we use the output BAM file from the previous process.  
- Then we sort the BAM file by read name using `samtools sort`. We will need this for the posterior processing.  

---

## Tools

| Tool     | Description                                   | Homepage |
|----------|-----------------------------------------------|----------|
| samtools | Alignment tool for mapping DNA sequences     | [GitHub](https://github.com/samtools/samtools) |

---

## Input

| Name   | Type   | Description          |
|--------|--------|--------------------|
| `sample` | string | Sample name         |
| `bam`    | file   | Path to the BAM file |

---

## Output

| Name       | Type | Description             |
|------------|------|-------------------------|
| `sample`     | string | Sample name          |
| `sorted_bam` | file   | BAM file sorted by read name |

---

# BAM_TO_ALLSITES_local

**Description:**  
Extracts potential viral integration sites from aligned BAM files.  
This module identifies putative integration sites by analyzing paired-end reads aligned against host genomic sequences. Stringent quality filters are applied, and integration coordinates are determined based on read orientation.  

The **start position** in the files corresponds to the R1 start in `+` strands and R2 start in `-` strands.

---

***Script details:***  

**1. Module and input loading**

- Loads the necessary modules and input files passed as arguments to the script.

**2. Quality filtering**

- Reads are filtered by proper pairing and mapping quality (`pysam` module):  
  [pysam AlignedSegment.is_proper_pair](https://pysam.readthedocs.io/en/stable/api.html#pysam.AlignedSegment.is_proper_pair)  
  Ensures INSPIIRED criteria:
  - Opposite strands
  - One read from each read (1 R1 + 1 R2)
  - Correct downstream orientation
- Removes duplicates, secondary, and supplementary alignments:  
  [is_duplicate](https://pysam.readthedocs.io/en/stable/api.html#pysam.AlignedSegment.is_duplicate)  
  [is_secondary](https://pysam.readthedocs.io/en/stable/api.html#pysam.AlignedSegment.is_secondary)  
  [is_supplementary](https://pysam.readthedocs.io/en/stable/api.html#pysam.AlignedSegment.is_supplementary)
- Filters by INSPIIRED criteria:
  - **Percent identity** threshold (calculated from CIGAR: `100 * matches / queryLength`)
  - **Soft-clipping limits** (sum soft-clipped bases at 5' end)
  - **Base insertion limits** (`tBaseInsert` ≤ 5bp)
  - **Alignment size** ≤ 2500 bp

    Creates a temporary file with all filtered individual reads (R1 and R2 separated).

**3. Paired-end processing**

- Pairs reads by line number and read name (BAM must be name-sorted).  
  - Unpaired reads or triplicates are discarded.
- Start position of insertions: R1 start on `+` strand, R2 start on `-` strand.  
  Captures the genomic breakpoint of viral integration.

**4. Output generation**

- Creates TSV with one entry per valid read pair.
- Includes `revmap` linking to original read line numbers in temporary file.
- Assigns sequential pairing IDs for tracking.

***Parameters:***
- `minPctIdent`: Ensures genuine alignments.
- `maxAlignStart`: Limits soft-clipping at 5' end.
- `tBaseInsert < 5`: Avoids large reference deletions.
- `maxFragLength`: Excludes chimeric fragments or mis-paired reads.

***Important notes:***
- BAM must be sorted by read name.
- Only properly paired reads contribute to integration sites.
- Each valid pair generates ONE integration site.
- Intermediate tmpFile contains filtered individual reads.
- Final `allSites.tsv` contains only successfully paired sites.

</details>

---

## Tools

<details>
<summary>python</summary>

**Description:** Python interpreter  
**Homepage:** [https://www.python.org/](https://www.python.org/)

</details>

<details>
<summary>pysam</summary>

**Description:** Python interface for SAM/BAM files  
**Homepage:** [https://pysam.readthedocs.io/](https://pysam.readthedocs.io/)

</details>

---

## Input

| Name | Type | Description | Pattern |
|------|------|-------------|---------|
| `sample` | string | Sample name | — |
| `bam` | file | Indexed and sorted BAM file | — |
| `minPctIdent` | string | Minimum percent identity (from SampleSheet) | — |
| `maxAlignStart` | string | Maximum soft-clipping at 5' end (from SampleSheet) | — |
| `maxFragLength` | string | Maximum fragment length (from SampleSheet) | — |

---

## Output

| Name | Type | Description | Pattern |
|------|------|-------------|---------|
| `sample` | string | Sample name | — |
| `*_allSites.tsv` | file | TSV with all insertion sites detected by INSPIIRED. Columns: `seqnames`, `start`, `end`, `strand`, `revmap`, `pairingID`, `samplename`, `ID` | `*_allSites.tsv` |
| `*_individualInsertions_notPaired.tsv` | file | TSV with all insertion sites not paired. Columns: `readname`, `chr`, `from`, `strand`, `start`, `end`, `cigar`, `qStart`, `PercIdent`, `tBaseInsert`, `flag`, `tags` | `*_individualInsertions_notPaired.tsv` |

---

# ALLSITES_TO_SITESFINAL_edited_grouping_local

**Description:**  
Insertions are initially deduplicated by collapsing all events that share the same **chromosome, start, end, and strand**. This ensures that identical insertions, which may arise from technical replication or multiple reads mapping to the same genomic coordinates, are represented only once. By doing so, we avoid artificially inflating counts and obtain a set of unique insertion events across the genome.  

The **counts** column in the metadata indicates the number of read pairs with identical insertion IDs, allowing assessment of duplicates per insertion. The **revmap** metadata contains indices of all reads mapping to the same insertion, enabling tracing back to the original reads.

**Function details:**  
The deduplication algorithm (`dereplicateSites_edited`) creates a unique identifier for each genomic position and collapses all reads mapping to that position. The `revmap` metadata allows tracing back to original reads for downstream analysis.

**Use case:**  
Accurately counts unique integration events and avoids inflating site numbers due to PCR duplicates or multiple reads from the same integration event.

---

## Tools

<details>
<summary>r-base</summary>

**Description:** R statistical computing environment  
**Homepage:** [https://www.r-project.org/](https://www.r-project.org/)

</details>

<details>
<summary>bioconductor-genomicranges</summary>

**Description:** Representation and manipulation of genomic intervals  
**Homepage:** [https://bioconductor.org/packages/GenomicRanges](https://bioconductor.org/packages/GenomicRanges)

</details>

<details>
<summary>bioconductor-hireadsprocessor</summary>

**Description:** Functions to process gene therapy vector integration sites  
**Homepage:** [https://bioconductor.org/packages/hiReadsProcessor](https://bioconductor.org/packages/hiReadsProcessor)

</details>

---

## Input

| Name | Type | Description | Pattern |
|------|------|-------------|---------|
| `sample` | string | Sample name | — |
| `allsites` | file | TSV file containing all insertion sites detected by INSPIIRED. Columns: `seqnames`, `start`, `end`, `strand`, `revmap`, `pairingID`, `samplename`, `ID` | — |
| `indivFile` | file | TSV file containing all insertion sites detected by INSPIIRED, not paired. Columns: `readname`, `chr`, `from`, `strand`, `start`, `end`, `cigar`, `qStart`, `PercIdent`, `tBaseInsert`, `flag`, `tags` | `*_individualInsertions_notPaired.tsv` |

---

## Output

| Name | Type | Description | Pattern |
|------|------|-------------|---------|
| `sample` | string | Sample name | — |
| `*_allsites.rds` | file | RDS file containing a GRanges object with all insertion sites detected by INSPIIRED. Columns: `seqnames`, `ranges`, `strand`, `revmap`, `pairingID`, `samplename` | `*_allsites.rds` |
| `*_sitesfinal.rds` | file | Deduplicated integration sites as GRanges object (PRIMARY output). Columns: `seqnames`, `ranges`, `strand`, `revmap`, `pairingID`, `samplename`, `counts` | `*_sitesfinal.rds` |

*Tip:* Use backticks (\``) to highlight file patterns or code elements.

---

## SITESFINAL_TO_POINTS_local

### Description

The insertions from allSites files are further collapsed based on their R2 start positions (viral insertion point), specifically, the *start for insertions on the + strand* and the *end for those on the - strand* (as the strand column always displays the read2 strand).  
This step collapses/joins insertions that occur at the same genomic point to quantify **how many times a same insertion occurs in different cells** (clonal expansion).  

This module adds a **new feature** to the existing INSPIIRED pipeline. As the original INSPIIRED pipeline was created for the analysis of the distribution of viral insertions, it is not possible to detect clonal expansion of the insertions.  
However, because the maximum insertion length is limited to 2500 bp (maxAlignmentLength), highly frequent insertions may be underestimated, as multiple insertions with identical R1 and R2 positions can happen by chance and be collapsed by the previous process.  
This means the resulting counts reflect site-specific occurrence, but insertions with a very high frequency should be interpreted as conservative estimates.

As an indicator of the clonal expansion level, we have the **counts** column, which is the number of reads supporting the same R2 start.  
Also, the **revmap** column allows us to identify the original insertions that form the clonal expansion.

**Script description:**

1. *Load libraries and arguments*  
   - Libraries: GenomicRanges, ChIPseeker, WriteXLS, rtracklayer, GenomicFeatures, GenomeInfoDb, dplyr, org.Hs.eg.db, clusterProfiler  
   - Arguments: sample, sitesFinal, reference genome, genome knownGene file  
   - Load the correct TxDb file and read the sitesFinal file

2. *Create output for empty files*  
   - If sitesFinal only contains a header line, create empty outputs for controls

3. *Collapse insertions*  
   - Convert GRanges object in sitesFinal RDS to dataframe  
   - Collapse insertions with the same R2 start positions:  
     - Create `r2_pos`:  
       - + strand → start position  
       - - strand → end position  
     - Group insertions by seqnames, strand, and r2_pos  
   - Create `counts` column (number of reads supporting the same R2 start)  
   - Create `revmap` column (original revmap information joined by `;`)

4. *Reconvert to GRanges*  
   - Convert collapsed dataframe back to GRanges object

5. *Annotate peaks*  
   - Annotate peaks with ChIPseeker `annotatePeak`  
   - tssRegion: -3000 to 3000  
   - TxDb: reference genome  
   - Annotation database: org.Hs.eg.db

6. *Generate plots*  
   - Coverage plot: `covplot`  
   - Pie chart: `plotAnnoPie`  
   - Bar chart: `plotAnnoBar`  
   - Upset plot: `upsetplot`  
   - Venn pie chart: `vennpie`  
   - Upset + venn pie chart: `upsetplot(vennpie=TRUE)`  
   - Distribution to TSS plot: `plotDistToTSS`

7. *Calculate GO and KEGG enrichment*  
   - Functions: `enrichGO`, `enrichKEGG` (clusterProfiler)  
   - Key type: ENTREZID  
   - Organism: hsa  
   - p-value cutoff: 0.05

8. *Generate Excel and PDF files*  
   - Excel: annotated peaks  
   - PDF: all plots (coverage, pie, bar, upset, venn pie, TSS distribution, GO and KEGG enrichment)

---

### Tools

| Tool | Description | Homepage |
|------|-------------|---------|
| GenomicRanges | Package for manipulating genomic intervals | [link](https://bioconductor.org/packages/release/bioc/html/GenomicRanges.html) |
| ChIPseeker | Package for analyzing ChIP-seq data | [link](https://bioconductor.org/packages/release/bioc/html/ChIPseeker.html) |
| WriteXLS | Package for writing Excel files | [link](https://cran.r-project.org/web/packages/WriteXLS/index.html) |
| rtracklayer | Package for reading and writing genome annotations | [link](https://bioconductor.org/packages/release/bioc/html/rtracklayer.html) |
| GenomicFeatures | Package for manipulating genomic features | [link](https://bioconductor.org/packages/release/bioc/html/GenomicFeatures.html) |
| GenomeInfoDb | Package for manipulating genome information | [link](https://bioconductor.org/packages/release/bioc/html/GenomeInfoDb.html) |
| dplyr | Package for manipulating data frames | [link](https://cran.r-project.org/web/packages/dplyr/index.html) |
| org.Hs.eg.db | Package for accessing gene annotation data | [link](https://bioconductor.org/packages/release/data/annotation/html/org.Hs.eg.db.html) |
| clusterProfiler | Package for clustering and annotating genes | [link](https://bioconductor.org/packages/release/bioc/html/clusterProfiler.html) |

---

### Input channel

| Name | Type | Description |
|------|------|-------------|
| `sample` | string | Sample name |
| `allSites` | file | Path to the allSites file |
| `sitesfinal` | file | Path to the sitesfinal file |
| `reference_genome_name` | string | Reference genome name |
| `reference genome` | file | Path to the reference genome file |
| `genome knownGene` | file | Path to the genome knownGene file |

---

### Output channel: points

| Name | Type | Description |
|------|------|-------------|
| `sample` | string | Sample name |
| `fig_points` | file | Path to the figure file |
| `annotated_points` | file | Path to the annotated file |

---
