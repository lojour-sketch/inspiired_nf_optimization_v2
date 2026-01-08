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

### Authors

- @liberentaizp
