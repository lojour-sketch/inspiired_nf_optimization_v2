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

## Authors

- @liberentaizp
