# pipeline_fixed — General changes vs pipeline

`pipeline_fixed` keeps the same overall INSPIIRED-derived workflow as `pipeline`, but introduces a small set of targeted changes to make the preprocessing entry points more robust and to expose additional QC information.

## 1. Preprocessing channel fix for BCL demultiplexing

The BCL preprocessing workflow now preserves the shape of the `BCL2FASTQ_local.out.fastq` output when passing reads to the UMI extraction setup.

- `ch_demux_fastq` is kept as the direct `(meta, [files])` tuple emitted by `BCL2FASTQ_local`
- Per-sample reads are regrouped explicitly before UMI extraction
- The auxiliary BCL outputs (`fastq_idx`, `undetermined`, `reports`, `stats`, `interop`) remain collected because they are only published or reused as grouped side outputs

This change prevents tuple-shape errors during downstream preprocessing while preserving the existing published outputs.

## 2. Unknown barcode QC after demultiplexing on both input paths

`pipeline_fixed` adds a dedicated `UNKNOWN_BARCODE_QC_local` step immediately after demultiplexing on both preprocessing branches: after `BCL2FASTQ_local` on the BCL path and after `DEMUXING_FASTQ_local` on the FASTQ path.

On the FASTQ path, the QC step uses assigned FASTQs after excluding `unmatched.*` outputs and uses the original `Undetermined_*_R*.fastq.gz` files from `--FASTQfolderDir` as the unknown-read input.

This step:
- estimates the fraction of reads assigned to samples versus undetermined reads
- reports the top unknown barcode sequences found in the undetermined R1 FASTQ
- classifies unknown barcodes as exact, near, or random relative to the normalized samplesheet
- summarizes low-diversity and G-homopolymer signals
- produces rule-based indicators to help distinguish barcode mismatch issues from random background or PhiX-like noise

The reports are published under `00_demux_unknown_qc/${params.projectName}`.

## 3. More robust annotation inside containers

`SITESFINAL_TO_POINTS_local` now handles `clusterProfiler` more safely in containerized environments.

- Cache-related environment variables are redirected to a writable task-local directory before KEGG enrichment runs
- KEGG enrichment is wrapped so that a failure in `enrichKEGG` no longer aborts the whole annotation process
- The remaining annotation outputs, plots, and tables are still produced even if KEGG enrichment is unavailable

## Summary

In short, `pipeline_fixed` differs from `pipeline` in three ways:
- safer BCL preprocessing channel wiring
- added demultiplexing unknown-barcode QC on both input paths
- more resilient annotation behavior in read-only container environments
