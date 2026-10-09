# Reference and environment configuration

`bushman.config` defines resources and reporting paths for the active original-compatible workflow. The default executor is Slurm; `-profile local` selects local execution. Override queues, resources and filesystem binds with a site-specific configuration such as `examples/cluster.config`.

Supply an explicit reference manifest based on `examples/references.json` or `reference_manifest.bushman.example.json`. The active alignment reference is a twoBit file, not a STAR index. Preserve the complete sequence set used by the original analysis, including alternative and patch sequences when present. Supply the frozen Bushman RefSeq annotation for annotation mode and, for complete original reference labels, the matching original Seqinfo RDS. Record file SHA256 values.

`container_manifest.json` identifies the two digest-pinned input-adapter containers. `containers/runtime.lock.json` identifies the separate frozen scientific SIF. The tested Nextflow version is 25.04.6. Source, effective metadata, runtime and reference identities participate in provenance and cache validation.

The legacy STAR reference auto-selection and resource overrides for removed processes are excluded from this release. They remain available in earlier Git commits.
