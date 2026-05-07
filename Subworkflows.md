# INSPIIRED subworkflows

# PREPROCESSING Subworkflow

This subworkflow performs preprocessing steps on sequencing data, including demultiplexing, UMI extraction, quality control, and filtering.

## Included Modules
- `NORMALIZE_index_length_local` from `../../../modules/local/normalize_index_length_local/main`
- If BCL2FASTQ demuxing: 
   - `BCL2FASTQ_local` from `../../../modules/local/bcl2fastq_local/main`
- If FASTQ demuxing: 
    - `DEMUXING_FASTQ_local` from `../../../modules/local/demuxing_fastq_local/main`
    - `CREATE_demux_samplesheet_local` from `../../../modules/local/create_demux_samplesheet_local/main`
- After either demultiplexing branch:
   - `UNKNOWN_BARCODE_QC_local` from `../../../modules/local/unknown_barcode_qc_local/main`
- `UMIEXTRACT_local` from `../../../modules/local/umi_extract_local/main`
- `FASTQCANDTRIM_wfl` from `../../../subworkflows/nf-core/fastqc_fastq_umitools_trimgalore_edited/main`
- `MULTIQC_wfl` from `../../../subworkflows/nf-core/multiqc/main`
- `LTRchecking_seqkit_local` from `../../../modules/local/LTR_primer_checking_local/main`
- `RCremoval_inspiired_local` from `../../../modules/local/rcremoval_local/main`
- `FINDVECTOR_local` from `../../../modules/local/find_vector_local/main`
- `SHORTREMOVE_local` from `../../../modules/local/short_seq_removal_local/main`

The description of the input and output channel structures of each module, along a more detailed description of the workflow, is provided in the specialized Markdown document for each module.

## Differences between BCL2FASTQ and FASTQ demuxing
The workflow differs depending on the input files:
    - If BCL2FASTQ demuxing:
        - The workflow uses the `BCL2FASTQ_local` module, which converts raw BCL files into per-sample FASTQ files.
      - Immediately after demultiplexing, the workflow runs `UNKNOWN_BARCODE_QC_local` on the assigned FASTQs and undetermined reads to summarize barcode quality and unknown-read composition.
    - If FASTQ demuxing:
        - The workflow does `CREATE_demux_samplesheet_local` to create a samplesheet with the barcodes that the posterior process will need.
            - The `--instrument` parameter is used only on this branch to decide whether `index2` must be reverse complemented before `fqtk` demultiplexing.
        - Then, the workflow does `DEMUXING_FASTQ_local` to demultiplex the FASTQ files.
      - Immediately after demultiplexing, the workflow runs `UNKNOWN_BARCODE_QC_local` on the assigned FASTQs excluding `unmatched.*` and on the original `Undetermined_*_R*.fastq.gz` inputs.

## Workflow: `PREPROCESSING_wfl`

### Inputs
| Channel | Description | Structure |
|---------|------|-------------|
| `ch_first_input` | Channel containing the first input tuple | Samplename, samplesheet, runfolder |
| `ch_linkers` | Channel containing sample linkers | Samplename, unique_linker, common_linker |
| `ch_primer_ltr` | Channel containing primer and LTR information | Samplename, primer, ltrbit, largeLTRFrag, mingDNA |
| `ch_vector` | Channel containing the vector sequence | vector fasta file |


### Main Steps

1. **Normalize Index Lengths**
   - Module: `NORMALIZE_index_length_local`
   - Inputs: `ch_first_input`
   - Outputs: `ch_bcl_input_normalized`, `ch_modified_samples`

   * Ensures all sample indexes have the same length, adjusting UMIs if necessary.

2. **Demultiplexing (BCL2FASTQ or FASTQ folder input)**
   - Modules: `BCL2FASTQ_local` or `DEMUXING_FASTQ_local`
   - Inputs: `ch_bcl_primer_input` or `ch_demux_input`
   - Outputs: `ch_demux_fastq`, `ch_fastq_idx`, `ch_undetermined`, `ch_undetermined_idx`, `ch_bcl_reports`, `ch_bcl_stats`, `ch_bcl_interop`

   * Converts raw BCL or undetermined FASTQ files into per-sample FASTQ files.

3. **Unknown barcode QC (BCL and FASTQ inputs)**
   - Module: `UNKNOWN_BARCODE_QC_local`
   - Inputs: assigned FASTQs from the active demultiplexing branch, unknown or undetermined FASTQs, normalized samplesheet
   - Published reports: `demux_unknown_barcode_qc.metrics.tsv`, `demux_unknown_barcode_qc.top_unknowns.tsv`, `demux_unknown_barcode_qc.indicators.txt`, `demux_unknown_barcode_qc.metrics.json`

   * Runs immediately after demultiplexing on both entry paths. On FASTQ input, the assigned FASTQ list excludes `unmatched.*` outputs and the unknown-read input comes from `--FASTQfolderDir/Undetermined_*_R*.fastq.gz`.
   * Summarizes the fraction of assigned versus unknown reads and classifies the unknown barcode population as exact, near, or random relative to the samplesheet.

4. **UMI Extraction**
   - Module: `UMIEXTRACT_local`
   - Inputs: `ch_umi_extract_input` (tuple of sample_id, linkers, reads, was_modified)
   - Outputs: `ch_umi_fastq`

   * Extracts UMIs from R1 reads and updates read headers.

5. **Quality Control and Trimming**
   - Subworkflow: `FASTQCANDTRIM_wfl`
   - Inputs: `ch_umi_out` (tuple of meta, reads)
   - Outputs: `fastqc_html`, `fastqc_zip`, `fastqc_html_trimmed`, `fastqc_zip_trimmed`, `ch_fastq_filtered`

   * Performs QC on reads and trims sequences as needed.

6. **MultiQC Reporting**
   - Subworkflow: `MULTIQC_wfl`
   - Inputs: `fastqc_html`, `fastqc_zip`, `fastqc_html_trimmed`, `fastqc_zip_trimmed`
   - Outputs: `multiqc_html`

   * Aggregates QC results from multiple samples.

7. **LTR and Primer Sequence Removal**
   - Module: `LTRchecking_seqkit_local`
   - Inputs: `ch_input` (tuple of sample_id, reads, primer, ltrbit, largeLTRFrag, project, mingDNA)
   - Outputs: `ch_ltr_chunks`

   * Removes reads that do not start with primer+LTR and cleans R2 sequences.

8. **Reverse Complement Removal**
   - Module: `RCremoval_inspiired_local`
   - Inputs: `ch_joined_input` (combined LTR chunks with linkers)
   - Outputs: `ch_rc_removed`

   * Removes reverse complement fragments from reads.

9. **Vector Sequence Filtering**
   - Module: `FINDVECTOR_local`
   - Inputs: `ch_findvector_input` (joined RC removed reads with primers) and `vectorfasta`
   - Outputs: `ch_vector_removed`

   * Identifies and removes vector sequences.

10. **Short Sequence Removal**
   - Module: `SHORTREMOVE_local`
   - Inputs: `ch_shortremove_input` (joined vector removed reads with primer info)
   - Outputs: `ch_short_removed`

   * Removes reads shorter than minimum length and their paired reads.



## Outputs
| Channel | Description | Structure |
|---------|------|-------------|
| `ch_short_removed` | Processed reads after all preprocessing steps | Samplename, R1 fastq file, R2 fastq file |

For both BCL and FASTQ inputs, the workflow also publishes demultiplexing QC reports to `${params.outdir}/00_demux_unknown_qc/${params.projectName}`. These reports are published side outputs rather than workflow channels.

## Notes
- The workflow handles two main entry points, and the selection of one is done depending on the --BCLorFASTQ parameter:
  1. Demuxing from BCL2FASTQ input.
  2. Demuxing from pre-existing FASTQ folder.
- Both paths run `UNKNOWN_BARCODE_QC_local` immediately after demultiplexing, then converge at UMI extraction and proceed through the same quality control and filtering steps.
- The FASTQ path uses `CREATE_demux_samplesheet_local`; this is the only preprocessing branch that consumes `--instrument`.
- All steps maintain paired-end read structure.
- The creation of new channels in between modules is done to alawys give a valid input channel to each process. As it is more robust to just give ONE channel as an input to each process, we join and manipulate the channels inside the workflow (and outside the modules) to ensure that the input channels are valid.

---

# ALIGNMENT Subworkflow

This subworkflow indexes the coresponding reference genome and aligns the reads to the reference genome.

---

## Included Modules

* `ALIGNMENT_local` from `'../../../modules/local/alignment/main'`
* `GENOME_INDEXING_local` from `'../../../modules/local/genome_indexing/main'`

The description of the input and output channel structures of each module, along a more detailed description of the workflow, is provided in the specialized Markdown document for each module.

---

## Workflow: `ALIGNMENT_wfl`

### Inputs
| Channel | Description | Structure |
|---------|-------------|-----------|
| `ch_short_removed` | Channel containing reads filtered for minimum length, output from preprocessing step. | Samplename, R1 fastq file, R2 fastq file |
| `ch_refGenome` | Channel containing reference genome information (name and paths from the samplesheet). | referenceGenome_name, referenceGenomeFile, referenceKnowngeneFile |

### Main Steps

1. **Get Unique Genomes for Indexing**

   ```groovy
   ch_unique_genomes = ch_refGenome
       .map { refGenome_name, refGenomeFile, refKnowngeneFile -> [refGenome_name, refGenomeFile] }
       .unique()
   ```

   * Extracts the genome name and reference fasta file.
   * Ensures only unique genomes are indexed, to not repeat an unnecessary indexing step.

2. **Run Genome Indexing**

   ```groovy
   GENOME_INDEXING_local(ch_unique_genomes)
   ```

   * Runs indexing for genomes that require it. If a genome is already indexed, it reuses the existing index, without running indexing again.

3. **Prepare Alignment Inputs**

   ```groovy
   ch_short_removed
       .merge(ch_refGenome)
       .map { sample, r1, r2, genome_name, refGenomeFile, refKnowngeneFile -> [genome_name, sample, r1, r2] }
       .set { ch_sample_with_genome }
   ```

   * Merges the filtered reads with reference genome information.

   ```groovy
   ch_sample_with_genome
       .combine(GENOME_INDEXING_local.out.index)
       .map { genome_name, sample, r1, r2, genome_name_duplicate, index_path -> [sample, r1, r2, index_path] }
       .set { ch_alignment_input }
   ```

   * Combines each sample with its genome index path for alignment.

5. **Run Alignment**

   ```groovy
   ALIGNMENT_local(ch_alignment_input)
   ALIGNMENT_local.out.aligned.set{ ch_aligned }
   ```

   * Performs alignment using the `ALIGNMENT_local` module.
   * Output channel `ch_aligned` contains the aligned reads.

### Outputs (Emit)

| Channel | Description | Structure |
|---------|-------------|-----------|
| `ch_aligned` | Channel containing aligned reads | Samplename, BAM file |

---

### Notes

* Multiple channels can be created for different types of alignments (chimeras, multimapped, uniquely mapped) if needed, as the alignment is run with STAR.
* This workflow ensures that each sample is aligned to its corresponding reference genome efficiently.
* Indexing is only performed for unique genomes to save computational resources.

---

# POSTPROCESSING Subworkflow

This workflow performs post-alignment processing including BAM sorting, generating allSites, collapsing insertions to sitesFinal, and creating points for clonal expansion analysis.

## Included Modules
- `INDEX_SORT_BAM_local` from `../../../modules/local/index_sort_bam_local/main`
- `BAM_TO_ALLSITES_local` from `../../../modules/local/bam_to_allsites_local/main` 
- `ALLSITES_TO_SITESFINAL_edited_grouping_local` from `../../../modules/local/allsites_to_sitesfinal_local/main`
- `SITESFINAL_TO_POINTS_local` from `../../../modules/local/sitesfinal_to_points_local/main`

The description of the input and output channel structures of each module, along a more detailed description of the workflow, is provided in the specialized Markdown document for each module.

## Workflow: `POSTPROCESSING_wfl`

### Inputs
| Channel | Description | Structure |
|---------|-------------|-----------|
| `ch_aligned` | Aligned BAM files from the alignment workflow | Samplename, BAM file |
| `ch_processing_params` | Parameters for post-processing (e.g., filtering thresholds) | minPctIdent, maxAlignStart, maxFragLength
| `ch_refGenome` | Reference genome information for annotation | referenceGenome_name, referenceGenomeFile, referenceKnowngeneFile |

### Main Steps
1. **BAM Sorting**
   - Process: `INDEX_SORT_BAM_local`
   - Input: `ch_aligned`
   - Output channel: `ch_sorted`

   * Sorts BAM files by coordinate and indexes them with samtools.

2. **Prepare BAM input**
   - Merge `ch_sorted` with `ch_processing_params`
   - Output channel: `ch_baminput`
   ```groovy
   ch_baminput = ch_sorted.merge(ch_processing_params)
   ```

3. **Convert BAM to allSites**
   - Process: `BAM_TO_ALLSITES_local`
   - Input: `ch_baminput`
   - Output channel: `ch_allsites`

   * Converts BAM files to allSites format, with one insertion per line.

4. **Collapse allSites to sitesFinal**
   - Process: `ALLSITES_TO_SITESFINAL_edited_grouping_local`
   - Input: `ch_allsites`
   - Output channel: `ch_sitesfinal`

    * Collapses allSites into sitesFinal by grouping duplicated insertions.

5. **Prepare annotation input**
   - Merge `ch_sitesfinal` with `ch_refGenome`
   - Output channel: `ch_annotinput`
   ```groovy
   ch_annotinput = ch_sitesfinal.merge(ch_refGenome)
   ```

6. **Collapse insertions and generate points**
   - Process: `SITESFINAL_TO_POINTS_local`
   - Input: `ch_annotinput`
   - Output channel: `ch_points`

   * Collapses insertions at the same genomic point (same R2 start, same viral insertion) to generate points for clonal expansion analysis. In `pipeline_fixed`, this step also redirects `clusterProfiler` cache usage to a task-local writable directory and continues if KEGG enrichment is unavailable.

### Outputs
| Channel | Description | Structure |
|---------|-------------|-----------|
| `ch_points` | Collapsed insertion points, counts, and clonal expansion information | Samplename, plots file, annotation excel

```
