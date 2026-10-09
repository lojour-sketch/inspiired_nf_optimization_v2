# Runtime and reproducibility

The executable environment is the exact tested `bushman_runtime.sif` image, distributed separately from the code. Its SHA256 is `4afa1dd705f7500c4b203a448cba06db35d580403270046e3a7a3a937418f4ed`; size is 1,697,886,208 bytes. The compact C routine compiles inside this image; gcc/make and the R ABI are included. Compiler output and native-source hashes are retained per caller task.

`containers/runtime.lock.json` identifies the image and tested headline versions. `vendor/bushman_r_sources/runtime_packages.tsv` records installed package versions. The larger runtime lock records constituent local-build files and their hashes. `containers/bushman_runtime.historical.def` is build-history evidence; its base image and separately compiled libraries are not in the small code package. The release therefore requires the tested SIF and does not claim a complete source-only rebuild.

The BCL/fixed FASTQ adapters use container image digests in nextflow.config. These images must be fetched or cached separately. Scientific references are external resources: record exact twoBit sequence content, SHA256, frozen RefSeq annotation, vector FASTA and optional original Seqinfo RDS. A new annotation database, genome subset, aligner or runtime is a new validation target.

To retain complete original caller reference metadata, export `seqinfo` from the original assembly-matched BSgenome or original saved mate-hit GRanges into an RDS, then provide `seqinfo_rds` and `seqinfo_sha256` in references.json. The caller validates exact twoBit sequence names, order and lengths before using those labels. Do not import a dictionary from another assembly or a subset reference.

Scientific source and effective metadata hashes participate in Nextflow task identity. Use fresh output directories for a corrected run; `-resume` is appropriate only with unchanged declared inputs or matching verified fingerprints. The old forward-marker RUN809 output is kept as historical evidence and does not become corrected output merely because the current source is fixed.
