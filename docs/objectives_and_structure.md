# Objective and relationship to the original Bushman pipeline

## Why develop this pipeline?

The objective is to make the original Bushman insertion-site analysis practical to install, run, resume, audit and maintain with current software. The scientific target remains the original analysis: identify vector–genome junctions, distinguish uniquely supported insertions from ambiguous mappings and chimeras, retain fragment evidence, estimate unique-site abundance with Sonic, and apply the original reporting and annotation functions.

The project is not a proposal for a new insertion-site caller. A faster implementation is useful only if it preserves the original evidence, candidate pairings, thresholds, classifications and downstream estimates. This is why the active workflow uses BLAT rather than adopting STAR as a speed shortcut, and why multihit processing and Sonic are retained. The optimized caller replaces expensive ways of locating and expanding rows; it does not intentionally replace the rules that decide what those rows mean biologically.

Nextflow addresses operational problems in the older deployment: manual stage submission, dependence on a particular server layout, completion-file bookkeeping, difficulty identifying the environment used for a result, and unnecessary waiting between samples. A frozen runtime addresses changing R, Bioconductor and system dependencies. Input and reference manifests make the effective analysis settings explicit. None of these improvements alone establishes that a result is biologically correct; they make the analysis reproducible and its scientific assumptions reviewable.

The release consequently has two acceptance criteria. First, the selected scientific outputs must agree with the committed original functions and the available historical original RUN809 evidence at a matched scope. Second, every deliberate difference must be identified as a compatibility adaptation, execution improvement, software repair, input extension or optional reporting choice. The release report records the completed checks and their limits rather than presenting all modernization changes as scientifically neutral by definition.

## The original analysis and its Nextflow counterpart

The following sequence describes the analysis that this release carries forward. The module-by-module report provides the exact files, changes and reasons for each operation.

```text
BCL or FASTQ input + sample metadata + frozen reference/vector inputs
    -> input conversion/routing and original index/Golay processing
    -> whole-sample read preparation and vector filtering
    -> sequence dereplication, paired-read keys and FASTA chunks
    -> independent BLAT alignment of both mates
    -> original paired-alignment tests and insertion-site classification
       -> unique insertion sites
       -> multihit connected components and their fragment evidence
       -> chimera classification and archive
    -> original report-group coordinate standardization
    -> original specimen grouping and fragment-based abundance
       -> Sonic for unique sites
       -> separate original length metric for multihit clusters
    -> frozen Bushman annotation and applicable population summaries
    -> exported evidence, status tables and additional review plots
```

| Analysis layer | Relationship to the original | Implementation in this release |
|---|---|---|
| Sample metadata and read assignment | Preserve the original effective metadata and index operations; provide explicit BCL/FASTQ entry points. | Input validation, BCL adapter and original Golay/index processing. The effective linker is derived from its complete sequence, including when a CSV supplies a conflicting marker. |
| Read preparation | Preserve the original whole-sample decisions, primer evidence, vector filtering, deduplication and paired-read keys. | A per-sample preparation process using frozen scientific statements and documented modern-header/leading-UMI adapters. |
| Alignment | Preserve independent mate alignments and the original BLAT parameters and reported hits. | Separate chunk tasks can run concurrently. A sample enters calling only when its complete expected alignment set is present. |
| Calling | Preserve candidate loci, admissible read-pair combinations, thresholds and unique/multihit/chimera classifications. | The default compact implementation computes memberships and intersections more efficiently. The original and indexed implementations remain selectable for comparison. |
| Ambiguous insertions | Preserve the original multihit connected-component logic and separate abundance evidence. | Ambiguous clusters are retained in their own objects, tables and plots; they are not assigned a guessed unique insertion coordinate. |
| Reporting coordinates | Preserve the original one-pass frequency/tie rules and report-group scope. | The original function has a configurable window; 5 bp is the default and 0 retains exact coordinates. Choosing another window is an analysis change. |
| Abundance | Preserve the original Sonic inputs, wrapper, rounding, proportions and ranks for unique sites. | Sonic is run per original specimen scope. Multihits retain their separate original length metric. Rejected Sonic input remains rejected. |
| Annotation and summaries | Carry forward the selected frozen Bushman functions and reference annotations. | Local frozen annotations replace dependence on a mutable installed reference. Applicable population calculations use the original functions with a recorded random seed. |
| Presentation | Preserve scientific outputs and expose them for review. | Added PDFs, annotated tables and an HTML index read existing output; they do not alter calling, standardization or estimation. |

The older R orchestration and the new Nextflow graph therefore differ structurally in how jobs are scheduled, inputs are staged, failures are tracked and completed work is resumed. Their intended scientific sequence is the same for the selected calling and abundance workflow. “Original-compatible” refers to that scientific contract, not to identical directory names, scripts, package versions or report appearance.

## What is intentionally different?

Current R/Bioconductor APIs, Python 3 syntax and modern Illumina headers require adapters. The RUN809 linker also begins with an N region, requiring the documented leading-UMI handling that was already used by the executable historical reference. These are compatibility changes; their exact scope is listed in the module table. They should not be confused with proof that every possible legacy library layout has been validated.

The chimera archive has an explicit bug repair. The historical loop counted data-frame columns rather than classified rows, saving four alignment entries for 14,760 chimera rows in the 0% archive. The corrected implementation saves all classified rows. The historical four entries match the corrected archive's prefix, and the insertion calls and classifications agree. This output is deliberately more complete, so a claim that every archive object is byte-for-byte identical would be incorrect.

Optional policies and analysis settings are exposed rather than hidden. Recording an unavailable Sonic estimate lets other valid specimens publish; it changes error handling, not the rejected sample's estimator input. Changing the site window changes reported coordinate grouping. Fragment-diversity abundance, fixed demultiplexing or disabling annotation are separate selectable modes, not evidence of original Sonic/index/annotation equivalence. The documented original-compatible controls should be used for comparison.

The new plotting module is a presentation addition. It uses frozen Bushman annotation and original abundance output to show insertion locations, annotation relationships and clonality evidence. It does not reproduce every ChIPseeker figure from the older Nextflow Step 17, nor does it reinterpret an unavailable estimate as zero.

## Scope within the wider INSPIIRED software suite

The original INSPIIRED project includes more than insertion calling. Its wider ecosystem also includes database upload and management, patient reporting, and genomic and epigenetic heatmap components. This release modernizes the selected insertion-calling, abundance, annotation and applicable population-report functions needed here. Direct RDS/TSV publication and the new review plots provide a transparent analysis record.

It does not recreate the complete SQL upload/database workflow, every original patient-report layout, longitudinal database management, genomic heatmaps or epigenetic heatmaps. Those are missing components relative to the wider suite, rather than removed scientific filters inside this caller. Reintroducing them would require their own input contracts, frozen dependencies and validation. In particular, this release should not be described as a complete replacement for the original clinical reporting system.

## What the current evidence establishes

The report includes a matched historical original 0% caller replay, matched single-specimen Sonic comparison, complete synthetic original/indexed/compact workflow comparisons, randomized operation oracles and final default-workflow checks. The status table records the separate full read-preparation check. Comparison of saved alignments tests the caller; it does not by itself test fresh raw-BCL conversion or a new genomic BLAT run.

The seven available earlier RUN809 outputs used an incorrect effective linker marker and remain explicitly provisional. They are useful for examining output structure and plots, but they are not the corrected release's seven-sample validation dataset. The corrected full seven-sample run, independent confirmed insertion coordinates, vector copy number and technical replicates remain necessary to assess mixture recovery and biological accuracy. The 0%/50%/75%/90%/95%/99%/100% labels describe percent polyclonal cells; an insertion's Sonic fraction is not automatically the fraction of cells in a clone containing multiple insertions.

This separation lets the repository explain both what has been preserved and what remains to be demonstrated. See the comprehensive HTML report for current evidence, the module table for precise modifications, and the results guide for interpreting the plots.
