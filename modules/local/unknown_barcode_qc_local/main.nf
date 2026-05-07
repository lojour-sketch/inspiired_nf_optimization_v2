process UNKNOWN_BARCODE_QC_local {

    publishDir "${params.outdir}/00_demux_unknown_qc/${params.projectName}", mode: 'copy', overwrite: true

    input:
    path(assigned_fastqs)
    path(unknown_fastqs)
    path(samplesheet)

    output:
    path("demux_unknown_barcode_qc.metrics.tsv"), emit: metrics
    path("demux_unknown_barcode_qc.top_unknowns.tsv"), emit: top_unknowns
    path("demux_unknown_barcode_qc.indicators.txt"), emit: indicators
    path("demux_unknown_barcode_qc.metrics.json"), emit: metrics_json

    script:
    """
    set -euo pipefail

    mkdir -p assigned_inputs unknown_inputs

    for f in ${assigned_fastqs}; do
        ln -s "\$(readlink -f "\$f")" assigned_inputs/ || true
    done

    for f in ${unknown_fastqs}; do
        ln -s "\$(readlink -f "\$f")" unknown_inputs/ || true
    done

    unknown_r1=\$(find unknown_inputs -maxdepth 1 \\( -type f -o -type l \\) -name "*R1*.f*q.gz" | head -n 1)
    if [[ -z "\$unknown_r1" ]]; then
        echo "ERROR: Could not find unknown/unmatched R1 FASTQ in unknown inputs" >&2
        ls -lah unknown_inputs >&2 || true
        exit 1
    fi

    mapfile -t assigned_files < <(find assigned_inputs -maxdepth 1 \\( -type f -o -type l \\) | sort)
    if [[ \${#assigned_files[@]} -eq 0 ]]; then
        echo "ERROR: No assigned FASTQ inputs found" >&2
        exit 1
    fi

    demux_unknown_barcode_qc.py \\
      --assigned-fastqs "\${assigned_files[@]}" \\
      --unknown-fastq "\$unknown_r1" \\
      --samplesheet ${samplesheet} \\
      --out-prefix demux_unknown_barcode_qc \\
      --top-n 10
    """
}
