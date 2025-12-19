# PREPROCESSING Subworkflow Documentation

This subworkflow performs preprocessing steps on sequencing data, including demultiplexing, UMI extraction, quality control, and filtering.

## Included Modules
- `NORMALIZE_index_length_local` from `../../../modules/local/normalize_index_length_local/main`
- If BCL2FASTQ demuxing: 
    - `BCL2FASTQ_local` from `../../../modules/local/bcl2fastq_local/main`
- If FASTQ demuxing: 
    - `DEMUXING_FASTQ_local` from `../../../modules/local/demuxing_fastq_local/main`
    - `CREATE_demux_samplesheet_local` from `../../../modules/local/create_demux_samplesheet_local/main`
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
    - If FASTQ demuxing:
        - The workflow does `CREATE_demux_samplesheet_local` to create a samplesheet with the barcodes that the posterior process will need.
        - Then, the workflow does `DEMUXING_FASTQ_local` to demultiplex the FASTQ files.

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

3. **UMI Extraction**
   - Module: `UMIEXTRACT_local`
   - Inputs: `ch_umi_extract_input` (tuple of sample_id, linkers, reads, was_modified)
   - Outputs: `ch_umi_fastq`

   * Extracts UMIs from R1 reads and updates read headers.

4. **Quality Control and Trimming**
   - Subworkflow: `FASTQCANDTRIM_wfl`
   - Inputs: `ch_umi_out` (tuple of meta, reads)
   - Outputs: `fastqc_html`, `fastqc_zip`, `fastqc_html_trimmed`, `fastqc_zip_trimmed`, `ch_fastq_filtered`

   * Performs QC on reads and trims sequences as needed.

5. **MultiQC Reporting**
   - Subworkflow: `MULTIQC_wfl`
   - Inputs: `fastqc_html`, `fastqc_zip`, `fastqc_html_trimmed`, `fastqc_zip_trimmed`
   - Outputs: `multiqc_html`

   * Aggregates QC results from multiple samples.

6. **LTR and Primer Sequence Removal**
   - Module: `LTRchecking_seqkit_local`
   - Inputs: `ch_input` (tuple of sample_id, reads, primer, ltrbit, largeLTRFrag, project, mingDNA)
   - Outputs: `ch_ltr_chunks`

   * Removes reads that do not start with primer+LTR and cleans R2 sequences.

7. **Reverse Complement Removal**
   - Module: `RCremoval_inspiired_local`
   - Inputs: `ch_joined_input` (combined LTR chunks with linkers)
   - Outputs: `ch_rc_removed`

   * Removes reverse complement fragments from reads.

8. **Vector Sequence Filtering**
   - Module: `FINDVECTOR_local`
   - Inputs: `ch_findvector_input` (joined RC removed reads with primers) and `vectorfasta`
   - Outputs: `ch_vector_removed`

   * Identifies and removes vector sequences.

9. **Short Sequence Removal**
   - Module: `SHORTREMOVE_local`
   - Inputs: `ch_shortremove_input` (joined vector removed reads with primer info)
   - Outputs: `ch_short_removed`

   * Removes reads shorter than minimum length and their paired reads.



## Outputs
| Channel | Description | Structure |
|---------|------|-------------|
| `ch_short_removed` | Processed reads after all preprocessing steps | Samplename, R1 fastq file, R2 fastq file |

## Notes
- The workflow handles two main entry points, and the selection of one is done depending on the --BCLorFASTQ parameter:
  1. Demuxing from BCL2FASTQ input.
  2. Demuxing from pre-existing FASTQ folder.
- Both paths converge at UMI extraction and proceed through the same quality control and filtering steps.
- All steps maintain paired-end read structure.
- The creation of new channels in between modules is done to alawys give a valid input channel to each process. As it is more robust to just give ONE channel as an input to each process, we join and manipulate the channels inside the workflow (and outside the modules) to ensure that the input channels are valid.