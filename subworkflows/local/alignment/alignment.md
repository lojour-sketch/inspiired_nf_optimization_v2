# ALIGNMENT Workflow

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
