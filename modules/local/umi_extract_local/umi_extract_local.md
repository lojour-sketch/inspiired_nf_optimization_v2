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

### Authors

- @liberentaizp
