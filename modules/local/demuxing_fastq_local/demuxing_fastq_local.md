# DEMUXING_FASTQ_local

## Description

Demultiplexes FASTQ files using `fastqk`.

This module identifies *undetermined reads* and demultiplexes them using a provided samplesheet and a specified read structure.

---

## Parameter Configuration

- **readStructure**  
  Defines the structure of the FASTQ reads and is passed as a command-line argument.  
  In INSPIIRED, the read structure is:

  `20B+T 12B +T`


    Meaning:
    - R1 barcode: first 20 bp of the read
    - R2 barcodes: located in the index files

- **sample-metadata**  
The samplesheet is provided to map each barcode to its corresponding sample. The SampleSheet format is described below.

---

## File Organization

- Demultiplexed FASTQ files are written to the corresponding sample folders under:

    'results/1_demuxed_undetermined/'


---

## Tools

| Tool | Description | Homepage |
|-----|------------|----------|
| fastqk | Demultiplexes FASTQ files | https://github.com/fulcrumgenomics/fqtk |

---

## Input

### Input overview

| Name | Type | Description | Pattern |
|-----|------|-------------|---------|
| `FASTQfolderDir` | path | Directory containing Undetermined FASTQ files | — |
| `samplesheet` | path | Illumina samplesheet in CSV format | `*.csv` |
| `readStructure` | string | FASTQ read structure | — |

---

### Samplesheet columns (detailed explanation)

The samplesheet must include the following columns:

| Column | Description |
|------|-------------|
| `Sample_ID` | Unique identifier for each sample |
| `index` | First index sequence (sample-unique linker in R1; first part of the iSL linker sequence, before the UMI) |
| `index2` | Second index sequence (Golay sequence in R2) |
| `common_linker` | Common linker for all R1 reads (last part of the iSL linker sequence, after the UMI) |
| `primer` | Primer sequence (R2) |
| `ltrbit` | LTR bit of the viral vector (R2) |
| `largeLTRFrag` | Large LTR fragment flag from the viral vector (R2) |
| `Sample_Project` | Project name used for data organization |
| `mingDNA` | Minimum genomic DNA threshold (30 bp in INSPIIRED) |
| `minPctIdent` | Minimum percent identity threshold (95 in INSPIIRED) |
| `maxAlignStart` | Maximum allowed soft-clipping at the 5′ end (5 bp in INSPIIRED) |
| `maxFragLength` | Maximum insertion fragment length, distance between R1 and R2 (2500 bp in INSPIIRED) |
| `refGenome` | Reference genome used for alignment (hg38 in INSPIIRED) |
| `vectorSeq` | Viral vector genome sequence path (`null` if insertion occurs via transposon) |

**Format notes:**
- The file must contain a `[Data]` section header
- The samplesheet path is provided via the command line

---

## Output

### Output overview

| Name | Type | Description | Pattern |
|-----|------|-------------|---------|
| `read1` | file | Demultiplexed FASTQ files for R1 reads | `*.R1.fq.gz` |
| `read2` | file | Demultiplexed FASTQ files for R2 reads | `*.R2.fq.gz` |

---

## Authors

- @liberentaizp
