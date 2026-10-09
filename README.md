# bushman_fast: original-compatible INSPIIRED in Nextflow

This release modernizes the original Bushman insertion-site analysis using Nextflow and a frozen runtime. It retains independent BLAT alignments, the original candidate-pair filters, multihit connected components, unique-site Sonic abundance, and the original coordinate-standardization rules. The default caller is the optimized `compact` implementation; `--caller_impl original` runs the frozen scientific expressions with the same compatibility adapters.

Read the [objective and original workflow structure](docs/objectives_and_structure.md), [module-by-module comparison](docs/pipeline_changes.md), [validation report](reports/pipeline_changes_report.html), and [results guide](docs/results.md). The report distinguishes compatibility adaptations, software bug fixes, execution improvements, and optional analysis changes. Historical RUN809 outputs created before the effective-linker correction are explicitly provisional and require rerunning.

## Install the tested environment

Use Linux x86_64, Nextflow 25.04.6 (the tested version), Java supported by that Nextflow version, and Singularity/Apptainer. Obtain the separate **bushman_runtime.sif** release artifact and verify it:

```bash
sha256sum /absolute/path/bushman_runtime.sif
# 4afa1dd705f7500c4b203a448cba06db35d580403270046e3a7a3a937418f4ed
```

The exact tested image contains R 4.4.1, sonicLength 1.4.7, and the scientific dependencies. The code archive does not contain the 1.6 GB image or genomic references. Until the image is uploaded separately, public users must obtain it from the release owner. The historical container recipe documents how the local image was assembled; it requires additional build inputs and is **not a complete portable rebuild recipe**. See [runtime and provenance](docs/runtime.md).

Prepare your samplesheet, FASTQ manifest, and reference manifest using [examples](examples/). Use absolute paths for vector FASTA, reference twoBit, frozen RefSeq annotation, optional original Seqinfo RDS, and FASTQ files. The twoBit must contain the same sequences as your original analysis; changing to a primary-chromosome-only reference changes multihit interpretation.

## Assigned FASTQs with original index processing

```bash
nextflow run main.nf -profile slurm \
  --BCLorFASTQ FASTQ --demux_mode original \
  --samplesheet /absolute/path/samples.csv \
  --fastq_manifest /absolute/path/fastqs.csv \
  --runfolderDir /absolute/path/run \
  --reference_manifest /absolute/path/references.json \
  --bushman_r_container /absolute/path/bushman_runtime.sif \
  --projectName MyRun --outdir /absolute/path/results/MyRun
```

The assigned FASTQ manifest must contain `Sample_ID,r1,index_fastq,r2`. The index FASTQ is required to reproduce original index-quality filtering and Golay correction. R1/R2 pairs alone are insufficient for this mode. The default index selection for RUN809 is I2; raw undetermined FASTQs use `--golay_index_read I1` or `I2` according to the actual library layout.

## BCL input

```bash
nextflow run main.nf -profile slurm \
  --BCLorFASTQ BCL --demux_mode original \
  --samplesheet /absolute/path/samples.csv \
  --runfolderDir /absolute/path/IlluminaRun \
  --reference_manifest /absolute/path/references.json \
  --bushman_r_container /absolute/path/bushman_runtime.sif \
  --bcl_bases_mask 'I20Y159,I12,Y143' --golay_index_read I2 \
  --projectName MyRun --outdir /absolute/path/results/MyRun
```

The BCL adapter uses a digest-pinned bcl2fastq container and preserves index reads for subsequent original processing. Configure the bases mask for your sequenced read layout; RUN809's mask is not universal. This adapter is inherited from the Nextflow modernization and is an input extension, rather than a reconstruction of every original BCL workflow. Its routing has synthetic validation; complete raw-BCL-to-historical equivalence is a separate validation scope.

## Important controls

| Option | Default | Meaning |
|---|---|---|
| `--caller_impl` | `compact` | Original scientific caller with faster row memberships/intersections; `original` and `indexed` remain available. |
| `--site_window_bp` | `5` | Original one-pass, frequency/tie standardization across a report group. `0` retains exact insertion coordinates. This does not change the multihit graph tolerance. |
| `--abundance_method` | `sonic` | Original unique-site estimator. `fragment_diversity` is explicitly a different mode. |
| `--abundance_failure_policy` | `fail` | Stop when Sonic rejects input. `record` publishes valid specimens and records missing abundance as unavailable without filtering invalid fragments. |
| `--annotation_mode` | `bushman` | Frozen Bushman RefSeq nearest-gene and overlap annotation; `off` disables it. |
| `--generate_plots` | `true` | Generate insertion-location, annotation, rank, cumulative-abundance and multihit PDFs with a browseable HTML index. |

Original/golay modes derive `linkerCommon` from `linkerSequence`, exactly as upstream metadata processing does. An explicitly supplied forward marker is overwritten and retained as `requested_linkerCommon` for auditing. For RUN809 the effective marker is **AGTCCCTTAAGCGGAG**. Optional fixed demultiplexing respects a supplied marker and is a separate processing mode.

The default configuration targets Slurm. Supply site-specific queues/resources with `-c your_cluster.config`. Use `-profile local` for a suitably sized workstation; RUN809-scale preparation/calling requests up to 96 GB. The synthetic validation uses smaller explicit resources.

## Test and inspect results

```bash
bash tests/run_validation.sh /absolute/path/bushman_runtime.sif
```

The suite checks effective metadata against the actual upstream R helper, randomized original-operation oracles, the original Sonic rejection rule, full synthetic original/indexed/compact caller and report parity, and the integrated plotting step. Full original RUN809 0% validation is recorded separately in the HTML report.

Open `results/MyRun/17_sitesfinal_to_points/MyRun/hg38__REPORT_GROUP/plots/index.html`. Scientific calls are in `04_original_calls`; modeled abundance, multihit abundance, estimator inputs and report status are in `05_original_report`. Read the [results guide](docs/results.md) before interpreting clonality.

Upstream scientific source revisions and hashes are recorded in `vendor/bushman/PROVENANCE.json`. License: GPL-3.0-or-later for the Bushman-derived code; third-party tools retain their own terms. See [third-party attribution](THIRD_PARTY.md).
