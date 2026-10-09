# bushman_fast 1.2.0 release notes

This release replaces the repository's active STAR-based modernization with an optimized Nextflow workflow whose scientific target is the original Bushman insertion-site analysis. Earlier implementation history remains in Git.

The active workflow retains independent original-parameter BLAT mate alignments, original candidate-pair filters, multihit connected components, unique-site Sonic abundance, separate multihit length evidence, frozen Bushman annotation and original report-group coordinate standardization. The default compact caller reduces repeated membership scans and row expansion. BCL and FASTQ input entry points remain available.

The effective linker metadata bug is corrected: original processing derives the reverse-complement suffix from the full linker sequence and overwrites a conflicting supplied marker. RUN809's effective marker is AGTCCCTTAAGCGGAG. All older seven-sample forward-marker results are labelled pre-correction and must be rerun for corrected biological comparison.

The corrected full historical 0% preparation matches all 85 FASTA chunks, read-pair keys, complete primer-ID contents and 12 preparation counters. Seven scientific caller payloads match the historical original, including complete reference labels when original Seqinfo is supplied. The deliberate chimera archive repair saves all 14,760 classified rows instead of four historical alignment entries; insertion calls and classification are unchanged. Matched single-specimen 5 bp Sonic reporting agrees at 74 sites and summed rounded abundance 99,038.

Synthetic workflows pass complete original/indexed/compact comparisons, randomized original-operation oracles and the Sonic input-domain failure test. The final default compact/local workflow completes with the new plotting process. Added clonality/location PDFs and HTML pages read existing outputs without changing scientific analysis.

The measured matched caller replay takes 555.174 seconds versus approximately 3 hours 7 minutes in the historical original. It excludes fresh genomic mapping and raw-input processing, and is not a full end-to-end benchmark. A complete corrected seven-sample RUN809 run and independent biological validation remain outstanding.

The exact tested scientific runtime is distributed separately as bushman_runtime.sif. The source package contains no raw reads, reference genomes, full RData archives or Nextflow work directory. It includes source, frozen upstream attribution/hashes, example inputs, tests, detailed explanations, compact validation evidence and labelled review reports. Obsolete active STAR modules and stale resource/reference configuration are excluded.

See objectives_and_structure.md, pipeline_changes.md, results.md and reports/pipeline_changes_report.html for the complete scope and interpretation.
