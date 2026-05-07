# pipeline_fixed — Changes vs pipeline

## Date
2026-04-30

## Background
The original pipeline (`pipeline/`) failed during a Run_606 execution. `PREPROCESSING_wfl:BCL2FASTQ_local` completed successfully (SLURM exit 0), but Nextflow immediately aborted with:

```
groovy.lang.MissingMethodException: No signature of method:
Script_d53e1443259175e7$..._closure17.call() is applicable for
argument types: (LinkedList) values: [[D131_d0_BCMA_BFP, ...R1_001.fastq.gz, ...]]
```

The crash happened in `subworkflows/local/preprocessing/main_bcl.nf`
(`Script_d53e1443259175e7` in the log), inside the `.map` operator that
feeds the UMI-extraction step.

---

## Root cause

In `main_bcl.nf`, after `BCL2FASTQ_local` finishes, its `.fastq` output
was collected with `.collect()`:

```groovy
// ORIGINAL (broken)
BCL2FASTQ_local.out.fastq.collect().set { ch_demux_fastq }
```

`BCL2FASTQ_local` is always run exactly once because its input is guarded
by `.take(1)`. This means `out.fastq` emits a **single** `(meta, [files])`
tuple — a 2-element structure.

`.collect()` gathers all items emitted by a channel into a single list.
Applied to a channel that emits one tuple, it wraps that tuple in another
list, producing a `LinkedList` like:

```
[[meta, [R1.fastq.gz, R2.fastq.gz]]]
```

The very next line tried to unpack it as a 2-element tuple:

```groovy
// This map expects (meta, files) — a 2-element tuple
.combine(ch_demux_fastq.map { meta, files -> files })
```

Groovy cannot call that closure with a single `LinkedList` argument, so
the workflow crashed before any downstream processes could be submitted.

---

## Fix applied

**File:** `subworkflows/local/preprocessing/main_bcl.nf`

Removed `.collect()` from the `ch_demux_fastq` assignment only. The other
outputs (`fastq_idx`, `undetermined`, `reports`, `stats`, `interop`) keep
`.collect()` because they are only stored/published and are never
destructured with a `map`.

```groovy
// FIXED
BCL2FASTQ_local.out.fastq                    .set { ch_demux_fastq }  // <-- .collect() removed
BCL2FASTQ_local.out.fastq_idx       .collect().set { ch_fastq_idx }
BCL2FASTQ_local.out.undetermined    .collect().set { ch_undetermined }
BCL2FASTQ_local.out.undetermined_idx.collect().set { ch_undetermined_idx }
BCL2FASTQ_local.out.reports         .collect().set { ch_bcl_reports }
BCL2FASTQ_local.out.stats           .collect().set { ch_bcl_stats }
BCL2FASTQ_local.out.interop         .collect().set { ch_bcl_interop }
```

With the fix, `ch_demux_fastq` holds the raw `(meta, [files])` tuple that
`BCL2FASTQ_local` emits, and the downstream `.map { meta, files -> files }`
destructures it correctly.

---

## How to run

Use `pipeline_fixed/main.nf` in place of `pipeline/main.nf`.
The run command (from `Run_606/nextflow_run606.sh`) becomes:

```bash
~/nextflow run /home/lojour/Bioinfo_projects/INSPIRED_chunhui/pipeline_fixed/main.nf \
  --BCLorFASTQ BCL \
  --genome hg38 \
  --samplesheet /home/lojour/Bioinfo_projects/INSPIRED_chunhui/Run_606/SampleSheet_run606.csv \
  --runfolderDir /home/lojour/Bioinfo_projects/INSPIRED_chunhui/Run_606/runfolder_606 \
  --outdir /home/lojour/Bioinfo_projects/INSPIRED_chunhui/Run_606/results_original \
  --projectName run606_Original \
  --readStructure '20B+T 12B +T' \
  --instrument "NextSeq2000" \
  -work-dir /home/lojour/Bioinfo_projects/INSPIRED_chunhui/Run_606/work_original \
  -resume
```

The `-resume` flag is safe to use: `NORMALIZE_index_length_local` and
`BCL2FASTQ_local` both completed successfully in the failed run and will
be pulled from the Nextflow cache.
