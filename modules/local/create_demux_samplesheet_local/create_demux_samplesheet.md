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

## Authors

- @lrenteria
