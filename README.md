# INSPIIRED optimization with NextFlow

This Nextflow pipeline is based in [INSPIIRED](https://github.com/BushmanLab/INSPIIRED?tab=readme-ov-file) software's [intSiteCaller](https://github.com/BushmanLab/intSiteCaller) module. It follows the same general workflow as the original pipeline, with several enhancements:
* Accepts BCL or FASTQ inputs
* Faster alignment and overall execution time
* Support for multiple samples
* An alternative insertion site detection strategy that counts each exact insertion once and enables identification of clonal expansions


## Workflow diagram

![Workflow_diagram](Workflow_image)


## Pre-requirements

In order to run this pipeline, some prerequisites must be met:
* `--samplesheet` parameter must be provided, which will be the path to the samplesheet. The samplesheet must contain the following columns:
  * `Sample_ID`: Sample ID
  * `index`: Index sequence (sample unique linker)
  * `ìndex2`: Second index sequence (Golay Sequence)
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
+ `--instrument` parameter must be provided, which will be the sequencing machine used.
* If the input is a BCL Run Folder:
  * `--runfolderDir` parameter must be provided, which will be the path to the BCL Run Folder. This folder must be in the same directory as the pipeline.
* If the input is a FASTQ file:
  * `--runfolderDir` parameter must be provided, which will be any folder in the same directory as the pipeline. It will be used as a reference path
  * `--FASTQfolderDir`parameter must be provided, which will be the path to the folder with Undetermined FASTQ files
  * `--readStructure`parameter must be provided, which will be the structure of template and barcode sequences. If a read has 34 nucleotides and the barcodes are of 12 nucleotides separately, the read structure is 34T 12B. However, if the barcodes are inside the read the read structure would be 12B34T. In our case, we will mostly have the following read structure: `20B+T 12B +T`
* The container images that are described in the .def files must be created and available.
* The fasta file of the vector's genomic sequence must be available.

## Running the pipeline

The pipeline can be run using the following command when running with a BCL input:

```
nextflow run main.nf \
    --BCLorFASTQ BCL \
    --runfolderDir /path/to/BCL/Run/Folder \
    --samplesheet /path/to/SampleSheet.csv \
    --projectName ProjectName \
    --readStructure '20B+T 12B +T' \
    -with-report reports/ProjectName_report.html \
    --instrument 'MiSeq' \
    -with-trace reports/ProjectName_trace.txt \
    -resume
```

When we have a FASTQ input, we can run the pipeline with the following command:

```
nextflow run main.nf \
    --BCLorFASTQ FASTQ \
    --runfolderDir /path/to/FASTQ/Run/Folder \
    --samplesheet /path/to/SampleSheet.csv \
    --FASTQfolderDir /path/to/Undetermined_FASTQ_Files \
    --readStructure '20B+T 12B +T' \
    --projectName ProjectName \
    -with-report reports/ProjectName_report.html \
    --instrument 'MiSeq' \
    -with-trace reports/ProjectName_trace.txt \
    -resume
```

## Output

The output of the pipeline will be in the `results` folder. The folder structure will be as follows:

```
results
├── 00_normalized_index_length
│   ├── ProjectName (every step will have a folder with the project name)
│   │   ├── files
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

## LICENSE

This project is licensed under the GPLv3 License - see the [LICENSE](LICENSE) file for details
