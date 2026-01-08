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

## Authors

- @lrenteria
