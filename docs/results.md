# Where results are and how to read them

The pipeline publishes a complete results directory under `--outdir`:

| Directory | Contents |
|---|---|
| `00_provenance/PROJECT/` | Validated effective sample metadata, requested/effective linker values, reference/runtime/source hashes. |
| `01_golay_demux/PROJECT/` | Reads selected by original index processing and its counters. |
| `02_original_trim/PROJECT/SAMPLE/` | Primer IDs, read-pair keys, trim statistics and curated mate FASTAs. |
| `03_original_blat/PROJECT/SAMPLE/` | Original-parameter BLAT PSL chunks and alignment timings. |
| `04_original_calls/PROJECT/SAMPLE/SAMPLE/` | Exact allSites, sites.final, multihit and chimera archives; caller metrics are one level above. |
| `05_original_report/PROJECT/ASSEMBLY__GROUP/` | Raw/standardized PCR fragment evidence, unique-site Sonic abundance/fit/input tables, separate multihit abundance, frozen annotation and estimator status. |
| `17_sitesfinal_to_points/PROJECT/ASSEMBLY__GROUP/plots/` | Browseable `index.html`, `fig_points_SPECIMEN.pdf`, `clonality_SPECIMEN.pdf`, `annotated_points_SPECIMEN.tsv`, status files and `clonality_summary.tsv`. |

Outputs are plotted by **specimen**, because Sonic uses the specimen's technical replicate strata. Multiple libraries belonging to one specimen are not presented as independent cell populations.

## Clonality figures

The abundance-rank plot shows dominant unique insertion sites. The cumulative plot shows how much estimated unique-site abundance is explained by the largest sites. Top-site bars and share plots identify coordinates to track across mixtures. Shannon diversity, Simpson concentration and effective-site count are descriptive statistics calculated from the already reported unique-site abundance; they do not replace the original populationInfo analyses.

The denominator is the sum of the original rounded **unique-site** abundance within the specimen. Multihits are excluded from that denominator and shown with their own original distinct-fragment-length metric. No mixture fraction is inferred directly from these plots. When supplied, `Expected_clonal_fraction` is biological context, not a target forced onto the estimator.

Location figures use strand-aware, 1-based insertion points and include the reference's alternative/fix sequences when present. Gene-overlap and nearest 5-prime feature figures use the existing frozen Bushman RefSeq annotations. The prior ChIPseeker exon/promoter pie, overlap UpSet/Venn and GO/KEGG enrichment pages are not reproduced using a different live database: those would be separate annotation analyses, not clonality measures or original-Bushman equivalence evidence.

When Sonic rejects a specimen, its report shows **unavailable**, retains insertion locations and multihit evidence, and does not show a fake abundance plot. The old pre-correction 95% sample fails the original 20 bp minimum because of a retained 2 bp fragment. Corrected-run behavior must be checked independently.

## Included RUN809 reports

`reports/RUN809_pre_correction/` contains plots and tables for all seven existing 8 October outputs, labelled provisional. They used the incorrect forward linker marker and must not be presented as corrected-release results. `reports/RUN809_corrected_0pct/`, when present, contains the corrected 0% original-alignment replay and its explicit validation scope. See the main HTML report for the latest completed checks.

## Matched original 0% comparison

The comparison uses the same complete archived original preparation and all 85 original genomic BLAT chunks. The report-group window is 5 bp and the abundance scope is the same single specimen. The original PCR evidence is recomputed through the frozen original reporting functions for this matched comparison.

| Quantity | Historical original evidence / matched original report | Corrected compact caller / matched report |
|---|---|---|
| Curated paired reads supplied to caller | 10,522,503 | 10,522,503 |
| Exact insertion sites before report standardization | 102 | 102 |
| Multihit clusters | 9 | 9 |
| Unique sites after original 5 bp standardization | 74 | 74 |
| Sum of rounded unique-site Sonic abundance | 99,038 | 99,038 |
| Largest unique-site Sonic fraction | Approximately 77.65% | Approximately 77.65% |

All seven scientific caller payloads pass complete object comparison, including hit reference labels when the original Seqinfo is supplied. The intentional exception is the repaired chimera archive: 14,760 classified rows are archived instead of the historical four alignment entries, with classification unchanged. Coordinate, abundance and rank fields agree exactly in the matched Sonic comparison; proportions agree within 1e-12 to allow TSV serialization.

The compact caller's measured worker time is 555.174 seconds, about 9 minutes 15 seconds. The historical original caller took about 3 hours 7 minutes. This comparison excludes new genomic mapping, input staging and complete raw-BCL processing; differences in hardware and R versions also apply. It supports the caller optimization, not a claim of a measured end-to-end speedup for a fresh seven-sample run.
