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

