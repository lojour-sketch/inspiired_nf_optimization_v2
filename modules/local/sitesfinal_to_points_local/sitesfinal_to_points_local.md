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

### Authors

- @liberentaizp
