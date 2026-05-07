# INSPIIRED optimization with Nextflow

This repository contains `pipeline_fixed`, a minimal Nextflow variant based on [liberentaizp/inspiired_nf_optimization](https://github.com/liberentaizp/inspiired_nf_optimization), which itself builds on [INSPIIRED](https://github.com/BushmanLab/INSPIIRED?tab=readme-ov-file) software's [intSiteCaller](https://github.com/BushmanLab/intSiteCaller) module.

Compared with the original Nextflow pipeline by liberentaizp, `pipeline_fixed` keeps the same overall workflow while adding some small fixes:
* Support for both BCL run-folder inputs and FASTQ-folder inputs
* A demultiplexing unknown-barcode QC report for both BCL and FASTQ runs
* More robust annotation in containers by using a writable local cache for `clusterProfiler` and continuing if KEGG enrichment is unavailable


## Prerequisites

In order to run this pipeline, some prerequisites must be met:
* The following top-level parameters are always required:
  * `--BCLorFASTQ`
  * `--samplesheet`
  * `--runfolderDir`
  * `--projectName`
  * `--outdir`
* `--samplesheet` must point to a sample sheet that contains the following columns:
  * `Sample_ID`: Sample ID
  * `index`: Index sequence (sample unique linker)
  * `index2`: Second index sequence (Golay Sequence)
  * `common_linker`: Common linker sequence
  * `primer`: Primer sequence
  * `ltrbit`: LTR bit sequence
  * `largeLTRFrag`: Large LTR fragment sequence
  * `Sample_Project`: Project name
  * `mingDNA`: Minimum DNA length
  * `minPctIdent`: Minimum percentage of identity
  * `maxAlignStart`: Maximum alignment start
  * `maxFragLength`: Maximum fragment length
  * `refGenome`: Reference genome name
  * `vectorSeq`: Vector sequence path
  
An example sample sheet is present in this repository as `Example_SampleSheet.csv`.

* If the input is a BCL Run Folder:
  * `--runfolderDir` must point to the BCL run folder.
  * `--instrument` is not required on this path. `bcl2fastq` handles index orientation internally.
  * `--readStructure` is not required on this path.
* If the input is a FASTQ folder:
  * `--runfolderDir` must still be provided, but it is only used as the run reference directory.
  * `--FASTQfolderDir` must point to the folder containing the undetermined FASTQ files.
  * `--readStructure` must describe the structure of template and barcode sequences. If a read has 34 template nucleotides and the barcodes are separate 12 nt reads, the read structure is `34T 12B`. If the barcodes are embedded in the read, the read structure could be `12B34T`. In this workflow the most common value is `20B+T 12B +T`.
  * `--instrument` only matters on this FASTQ path, because `CREATE_demux_samplesheet_local` uses it to decide whether `index2` stays forward-oriented or is reverse-complemented before `fqtk` demultiplexing. Supported values are `MiSeq`, `NextSeq2000`, and `NextSeq500`.
* The container images described in the `.def` files must be created and available.
* The FASTA file of the vector's genomic sequence must be available in the same directory as the pipeline.
* The FASTA file of the reference genome must be available in the same directory as the pipeline, and its name must start with the genome name (hg19, hg38...) and finish with the `.fa` extension.

## Running the pipeline

The pipeline can be run using the following command when running with a BCL input:

```
nextflow run main.nf \
    --BCLorFASTQ BCL \
    --runfolderDir /path/to/BCL/Run/Folder \
    --samplesheet /path/to/Example_SampleSheet.csv \
    --projectName ProjectName \
    --outdir /path/to/results \
    -with-report reports/ProjectName_report.html \
    -with-trace reports/ProjectName_trace.txt \
    -resume
```

When using a FASTQ input folder, we can run the pipeline with the following command:

```
nextflow run main.nf \
    --BCLorFASTQ FASTQ \
    --runfolderDir /path/to/FASTQ/Run/Folder \
    --samplesheet /path/to/Example_SampleSheet.csv \
    --FASTQfolderDir /path/to/Undetermined_FASTQ_Files \
    --readStructure '20B+T 12B +T' \
    --projectName ProjectName \
    --outdir /path/to/results \
    -with-report reports/ProjectName_report.html \
    --instrument 'NextSeq2000' or 'MiSeq' \
    -with-trace reports/ProjectName_trace.txt \
    -resume
```

## Output

The output of the pipeline is written under `--outdir`. A typical folder structure is:

```
results
├── 00_create_demux_samplesheet
│   └── ProjectName                  # FASTQ input only
├── 00_demux_unknown_qc
│   └── ProjectName                  # BCL and FASTQ inputs
├── 00_normalized_index_length
│   └── ProjectName                  # every published step uses a project-specific subfolder
├── 1_demuxed
├── 2_extractedumi
├── 3_fastqcraw
├── 4_trimmedfastq_fastqc
├── 5_removed_n
├── 6_fastqctrimmed
├── 7_multiqcaftertrim
├── 8_LTR_presence
├── 9_reverse_complement_removal
├── 10_findvector
├── 11_short_remove
├── 12_genome_index
├── 13_alignment
├── 14_index_sort_bam
├── 15_allsites
├── 16_sitesfinal
├── 17_sitesfinal_to_points

```

For both BCL and FASTQ inputs, `00_demux_unknown_qc/ProjectName` contains four summary files:
* `demux_unknown_barcode_qc.metrics.tsv`
* `demux_unknown_barcode_qc.top_unknowns.tsv`
* `demux_unknown_barcode_qc.indicators.txt`
* `demux_unknown_barcode_qc.metrics.json`

On the FASTQ path, this QC step uses the demultiplexed sample FASTQs together with the original `Undetermined_*` FASTQs from `--FASTQfolderDir`.

## LICENSE

This project is licensed under the GPLv3 License - see the [LICENSE](LICENSE) file for details
