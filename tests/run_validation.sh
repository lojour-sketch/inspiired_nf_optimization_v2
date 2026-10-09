#!/usr/bin/env bash
set -euo pipefail
release_root=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)
runtime=${1:?Usage: tests/run_validation.sh /absolute/path/bushman_runtime.sif [validation_directory]}
test_root=${2:-"$release_root/results/synthetic_validation"}
mkdir -p "$test_root" "$test_root/fixture"
test_root=$(cd -- "$test_root" && pwd)
runtime=$(realpath "$runtime")
export NXF_HOME="$test_root/nextflow_home"
export NXF_OPTS='-Xms256m -Xmx3g'
export NXF_SINGULARITY_CACHEDIR="$test_root/container_cache"
mkdir -p "$NXF_HOME" "$NXF_SINGULARITY_CACHEDIR"
cd "$test_root"
run_r() { singularity exec --cleanenv --bind "$release_root,$test_root" "$runtime" Rscript "$@"; }
python3 "$release_root/tests/test_effective_linker.py" "$release_root" "$test_root/linker_cases.json"
run_r "$release_root/tests/test_effective_linker.R" "$release_root" "$test_root/linker_cases.json" "$test_root/effective_linker_parity.json"
run_r "$release_root/tests/test_indexed_kernels.R" "$release_root" "$test_root"
run_r "$release_root/tests/test_compact_pairing.R" "$release_root" "$test_root"
run_r "$release_root/tests/test_abundance_failure_policy.R" "$release_root/bin/bushman_pool.R" "$test_root/abundance_policy.json"
run_r "$release_root/tests/make_bushman_fixtures.R" "$test_root/fixture" "$release_root"
singularity exec --cleanenv --bind "$release_root,$test_root" "$runtime" python3 \
    "$release_root/tests/make_shared_core_fixture.py" "$test_root/fixture" "$release_root"
cat > "$test_root/local.config" <<CFG
process {
    executor='local';queue=null;cpus=1;memory='4GB';maxForks=2
    withName:BUSHMAN_GOLAY_SAMPLE { memory='4GB';cpus=1 }
    withName:BUSHMAN_TRIM { memory='4GB';cpus=1 }
    withName:BUSHMAN_BLAT { memory='4GB';cpus=1;maxForks=2 }
    withName:BUSHMAN_CALL { memory='4GB';cpus=1 }
    withName:BUSHMAN_POOL { memory='4GB';cpus=1 }
    withName:BUSHMAN_PLOTS { memory='4GB';cpus=1 }
}
singularity.runOptions='--cleanenv --bind $release_root,$test_root'
CFG
for implementation in original indexed compact; do
    nextflow -log "$test_root/$implementation.nextflow.log" run "$release_root/main.nf" \
        -c "$test_root/local.config" -work-dir "$test_root/${implementation}_work" -ansi-log false \
        --BCLorFASTQ FASTQ --demux_mode original --golay_index_read I2 --caller_impl "$implementation" \
        --samplesheet "$test_root/fixture/samples.csv" --runfolderDir "$test_root/fixture" \
        --fastq_manifest "$test_root/fixture/fastqs.csv" --reference_manifest "$test_root/fixture/references.json" \
        --bushman_r_container "$runtime" --annotation_mode off --projectName truth \
        --outdir "$test_root/${implementation}_results" --abundance_method sonic --site_window_bp 5
done
for implementation in indexed compact; do
    run_r "$release_root/tests/compare_payloads.R" "$test_root/original_results" truth \
        "$test_root/${implementation}_results" truth "$test_root/${implementation}_payload_parity.json"
done
echo 'PASS: effective metadata, randomized pairing, complete caller/report parity and Nextflow plotting'
