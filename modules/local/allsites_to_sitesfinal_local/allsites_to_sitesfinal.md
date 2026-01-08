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

## Authors

- @liberentaizp
