def bushmanImplementationHash(root, stage) {
    if(stage=='source') {
        def files = []
        ['bin','modules','subworkflows','conf','vendor/bushman'].each { folder ->
            new File("${root}/${folder}").eachFileRecurse { f ->
                if(f.isFile() && !f.path.contains('__pycache__')) files.add(f)
            }
        }
        ['main.nf','nextflow.config','pipeline_identity.json','vendor/bushman_r_sources/runtime_lock.json'].each { name -> files.add(new File("${root}/${name}")) }
        def sourceDigest = java.security.MessageDigest.getInstance('SHA-256')
        files.sort { it.path }.each { f -> sourceDigest.update(f.path.substring(root.toString().length()).getBytes('UTF-8'));sourceDigest.update(f.bytes) }
        return sourceDigest.digest().encodeHex().toString()
    }
    def dependencies = [
        validate: ['bushman_inputs.py'],
        trim: ['bushman_trim.R', 'bushman_runtime.R'],
        golay: ['bushman_golay.R', 'bushman_golay.py', 'bushman_runtime.R'],
        call: ['bushman_call.R', 'bushman_runtime.R', 'bushman_indexed.R', 'bushman_compact.R', 'bushman_pairing.c'],
        pool: ['bushman_pool.R', 'bushman_runtime.R'],
        plots: ['bushman_plots.R']
    ]
    if(!dependencies.containsKey(stage)) throw new IllegalArgumentException("Unknown Bushman stage: ${stage}")
    def digest = java.security.MessageDigest.getInstance('SHA-256')
    dependencies[stage].sort().each { name ->
        digest.update(name.getBytes('UTF-8'))
        digest.update(new File("${root}/bin/${name}").bytes)
    }
    def lock = new File("${root}/vendor/bushman_r_sources/runtime_lock.json")
    if(lock.exists()) digest.update(lock.bytes)
    digest.update(new File("${root}/vendor/bushman/PROVENANCE.json").bytes)
    digest.update(new File("${root}/pipeline_identity.json").bytes)
    digest.digest().encodeHex().toString()
}

process BUSHMAN_SOURCE_PROVENANCE {
    publishDir "${params.outdir}/00_provenance/${params.projectName}", mode: 'copy'
    input:
    val implementation
    output:
    path 'pipeline_implementation.json', emit: provenance
    script:
    """
    export BUSHMAN_IMPLEMENTATION_SHA256='${implementation}'
    python3 '${projectDir}/bin/bushman_source_manifest.py' '${projectDir}' \
        ${params.caller_impl} ${params.site_window_bp} ${params.abundance_method} ${params.annotation_mode}
    """
}

process BUSHMAN_VALIDATE {
    publishDir "${params.outdir}/00_provenance/${params.projectName}", mode: 'copy'
    input:
    path sheet
    path references
    path fastq_manifest
    path sources
    val implementation
    output:
    path 'validated_samples.tsv', emit: metadata
    path 'samples/*.json', emit: rows
    path 'input_provenance.json', emit: provenance
    script:
    def manifestArg = params.fastq_manifest ? "--fastq-manifest ${fastq_manifest}" : ''
    """
    export BUSHMAN_IMPLEMENTATION_SHA256='${implementation}'
    bushman_inputs.py ${sheet} '${params.runfolderDir}' ${references} '${projectDir}' \
        ${manifestArg} --annotation ${params.annotation_mode} --demux ${params.demux_mode} --input-type ${params.BCLorFASTQ} \
        --runtime-image '${params.bushman_r_container}' --blat-binary '${params.blat_binary}'
    """
}

process BUSHMAN_TRIM {
    tag "$sample"
    publishDir "${params.outdir}/02_original_trim/${params.projectName}/${sample}", mode: 'copy'
    input:
    tuple val(sample), path(row), path(r1), path(r2), path(vector), path(twobit), val(genome)
    path sources
    path(blat, stageAs: 'blat')
    val implementation
    output:
    tuple val(sample), path(row), path("${sample}"), path(twobit), val(genome), emit: prepared
    tuple val(sample), path("${sample}/R[12]-*.fa"), path(twobit), val(genome), optional: true, emit: chunks
    tuple val(sample), path('trim_metrics.json'), path('trim_runtime.json'), emit: audit
    script:
    """
    export BUSHMAN_IMPLEMENTATION_SHA256='${implementation}'
    Rscript '${projectDir}/bin/bushman_trim.R' ${row} ${r1} ${r2} ${vector} ${sources} ./blat ${params.blat_chunk_size}
    """
}

process BUSHMAN_GOLAY {
    publishDir "${params.outdir}/01_golay_demux/${params.projectName}", mode: 'copy'
    input:
    path metadata
    path folder
    path sources
    val implementation
    output:
    path '*.R[12].fq.gz', emit: reads
    path 'golay_runtime.json', emit: audit
    path 'golay_metrics.json', emit: metrics
    script:
    """
    export BUSHMAN_IMPLEMENTATION_SHA256='${implementation}'
    set -euo pipefail
    mapfile -t r1 < <(find -L ${folder} -maxdepth 1 -name 'Undetermined_*_R1*.fastq.gz' -type f | sort)
    mapfile -t ix < <(find -L ${folder} -maxdepth 1 -name 'Undetermined_*_${params.golay_index_read}*.fastq.gz' -type f | sort)
    mapfile -t r2 < <(find -L ${folder} -maxdepth 1 -name 'Undetermined_*_R2*.fastq.gz' -type f | sort)
    [[ \${#r1[@]} -eq 1 && \${#ix[@]} -eq 1 && \${#r2[@]} -eq 1 ]]
    Rscript '${projectDir}/bin/bushman_golay.R' ${metadata} "\${r1[0]}" "\${ix[0]}" "\${r2[0]}" ${sources}
    """
}

process BUSHMAN_GOLAY_SAMPLE {
    tag "$sample"
    publishDir "${params.outdir}/01_golay_demux/${params.projectName}/${sample}", mode: 'copy'
    input:
    tuple val(sample), path(row), path(r1, stageAs: 'input_R1.fastq.gz'), path(index_read, stageAs: 'input_index.fastq.gz'), path(r2, stageAs: 'input_R2.fastq.gz')
    path sources
    val implementation
    output:
    tuple val(sample), path("${sample}.R1.fq.gz"), path("${sample}.R2.fq.gz"), emit: reads
    tuple val(sample), path('golay_metrics.json'), path('golay_runtime.json'), emit: audit
    script:
    """
    export BUSHMAN_IMPLEMENTATION_SHA256='${implementation}'
    Rscript '${projectDir}/bin/bushman_golay.R' ${row} ${r1} ${index_read} ${r2} ${sources}
    """
}

process BUSHMAN_BLAT {
    tag "${sample}:${query.name}"
    publishDir "${params.outdir}/03_original_blat/${params.projectName}/${sample}", mode: 'copy'
    input:
    tuple val(sample), path(query), path(twobit), val(genome), val(completion_key)
    path(blat, stageAs: 'blat')
    output:
    tuple val(completion_key), path("${query.name}.psl.gz"), emit: psl
    tuple val(sample), path("${query.name}.timing.txt"), emit: timing
    script:
    """
    { time -p ./blat ${twobit} ${query} ${query.name}.psl \
        -tileSize=11 -stepSize=9 -minIdentity=85 -maxIntron=5 -minScore=27 -dots=1000 -out=psl -noHead; } 2> ${query.name}.timing.txt
    gzip ${query.name}.psl
    """
}

process BUSHMAN_CALL {
    tag "$sample"
    publishDir "${params.outdir}/04_original_calls/${params.projectName}/${sample}", mode: 'copy'
    input:
    tuple val(sample), path(row), path(prepared, stageAs: 'prepared'), path(twobit), val(genome), path(psls, stageAs: 'psl/*')
    path sources
    val implementation
    output:
    tuple val(sample), path("${sample}"), val(genome), path(row), emit: calls
    tuple val(sample), path('calling_metrics.json'), path('calling_runtime.json'), emit: audit
    script:
    """
    export BUSHMAN_IMPLEMENTATION_SHA256='${implementation}'
    mkdir -p psl
    Rscript '${projectDir}/bin/bushman_call.R' ${row} ${prepared} ${twobit} ${sources} psl ${params.caller_impl}
    """
}

process BUSHMAN_POOL {
    tag "$genome"
    publishDir "${params.outdir}/05_original_report/${params.projectName}/${genome}", mode: 'copy'
    input:
    tuple val(genome), path(calls, stageAs: 'calls/*')
    path metadata
    path sources
    val implementation
    output:
    tuple val(genome), path('original_*'), path('pooled_sites.rds'), path('sample_stats.tsv'), path('abundance_status.tsv'), emit: results
    tuple val(genome), path('pooled_metrics.json'), path('pooled_runtime.json'), emit: audit
    script:
    """
    export BUSHMAN_IMPLEMENTATION_SHA256='${implementation}'
    Rscript '${projectDir}/bin/bushman_pool.R' ${metadata} calls ${sources} \
        ${params.site_window_bp} ${params.abundance_method} ${params.annotation_mode} ${genome} ${params.abundance_failure_policy}
    """
}

process BUSHMAN_PLOTS {
    tag "$genome"
    publishDir "${params.outdir}/17_sitesfinal_to_points/${params.projectName}/${genome}", mode: 'copy'
    input:
    tuple val(genome), path(report_files, stageAs: 'report/*'), path(pooled_sites, stageAs: 'report/pooled_sites.rds'), path(sample_stats, stageAs: 'report/sample_stats.tsv'), path(abundance_status, stageAs: 'report/abundance_status.tsv')
    path metadata
    val implementation
    output:
    tuple val(genome), path('plots'), emit: figures
    script:
    """
    export BUSHMAN_IMPLEMENTATION_SHA256='${implementation}'
    Rscript '${projectDir}/bin/bushman_plots.R' report ${metadata} ${params.site_window_bp} plots '${params.projectName}'
    """
}
