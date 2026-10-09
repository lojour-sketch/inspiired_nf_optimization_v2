include { BUSHMAN_VALIDATE; BUSHMAN_TRIM; BUSHMAN_BLAT; BUSHMAN_CALL; BUSHMAN_POOL; BUSHMAN_PLOTS; BUSHMAN_GOLAY; BUSHMAN_GOLAY_SAMPLE; BUSHMAN_SOURCE_PROVENANCE; bushmanImplementationHash } from '../../../modules/local/bushman/main'
include { NORMALIZE_index_length_local } from '../../../modules/local/normalize_index_length_local/main'
include { CREATE_demux_samplesheet_local } from '../../../modules/local/create_demux_samplesheet_local/main'
include { DEMUXING_FASTQ_local } from '../../../modules/local/demuxing_fastq_local/main'
include { BCL2FASTQ_local } from '../../../modules/local/bcl2fastq_local/main'

workflow ORIGINAL_COMPATIBLE {
    main:
    sources = Channel.value(file("${projectDir}/vendor/bushman",checkIfExists:true))
    blat = Channel.value(file(params.blat_binary,checkIfExists:true))
    manifest = Channel.value(file(params.fastq_manifest ?: "${projectDir}/conf/empty_fastq_manifest.csv",checkIfExists:true))
    // Hashes are declared value inputs, computed before task cache lookup.
    // Computing them inside a process script closure is too late for caching.
    stageHashes = ['validate','trim','golay','call','pool','plots','source'].collectEntries { stage ->
        [(stage):Channel.value(bushmanImplementationHash(projectDir,stage))]
    }
    BUSHMAN_SOURCE_PROVENANCE(stageHashes.source)
    BUSHMAN_VALIDATE(Channel.of(file(params.samplesheet,checkIfExists:true)),
        Channel.value(file(params.reference_manifest,checkIfExists:true)),manifest,sources,stageHashes.validate)
    metadata = BUSHMAN_VALIDATE.out.metadata.first()
    rows = BUSHMAN_VALIDATE.out.rows.flatten().map { json ->
        def row = new LinkedHashMap(new groovy.json.JsonSlurper().parseText(json.text))
        tuple(row.Sample_ID, json, row)
    }
    if (params.fastq_manifest) {
        if(params.demux_mode=='original') {
            golayInput = rows.map { sample, json, row -> tuple(sample,json,
                file(row.r1,checkIfExists:true),file(row.index_fastq,checkIfExists:true),file(row.r2,checkIfExists:true)) }
            BUSHMAN_GOLAY_SAMPLE(golayInput,sources,stageHashes.golay)
            ingress = BUSHMAN_GOLAY_SAMPLE.out.reads
        } else {
            ingress = rows.map { sample, json, row -> tuple(sample,file(row.r1,checkIfExists:true),file(row.r2,checkIfExists:true)) }
        }
    } else if (params.demux_mode=='golay' || (params.demux_mode=='original' && params.BCLorFASTQ=='FASTQ')) {
        BUSHMAN_GOLAY(metadata,Channel.value(file(params.FASTQfolderDir,type:'dir',checkIfExists:true)),sources,stageHashes.golay)
        ingress = BUSHMAN_GOLAY.out.reads.flatten()
            .map { read -> tuple(read.name.replaceFirst(/\.R[12]\.fq\.gz$/, ''),read) }
            .groupTuple(by:0).map { sample, pair ->
                tuple(sample,pair.find { it.name.endsWith('.R1.fq.gz') },pair.find { it.name.endsWith('.R2.fq.gz') })
            }
    } else {
        first = rows.first().map { sample, json, row ->
            tuple([id:'run1'],file(params.samplesheet,checkIfExists:true),file(params.runfolderDir,type:'dir',checkIfExists:true))
        }
        NORMALIZE_index_length_local(first)
        if (params.BCLorFASTQ=='BCL') {
            bclInput = NORMALIZE_index_length_local.out.normalized.combine(rows.first()).map { meta, sheet, folder, sample, json, row ->
                tuple(sample,row.primer,row.ltrbit,row.largeLTRFrag,row.Sample_Project,row.mingDNA,meta,sheet,folder)
            }
            BCL2FASTQ_local(bclInput)
            bclFiles = params.demux_mode=='original' ?
                BCL2FASTQ_local.out.fastq.join(BCL2FASTQ_local.out.fastq_idx,by:0,failOnDuplicate:true,failOnMismatch:true)
                    .map { meta,reads,index_reads -> tuple(meta,
                        (reads instanceof List ? reads : [reads])+(index_reads instanceof List ? index_reads : [index_reads])) } :
                BCL2FASTQ_local.out.fastq
            assigned = bclFiles.flatMap { meta, files ->
                def paired = files.findAll { it.name.contains('_R1_') || it.name.contains('_R2_') ||
                    (params.demux_mode=='original' && it.name.contains("_${params.golay_index_read}_")) }
                    .groupBy { it.name.replaceFirst(/_S[0-9]+_.*$/, '') }
                paired.collect { sample, pair ->
                    def count = params.demux_mode=='original' ? 3 : 2
                    if (pair.size()!=count) error "Expected synchronized reads and index for ${sample}"
                    if(params.demux_mode=='original')
                        tuple(sample,pair.find { it.name.contains('_R1_') },
                            pair.find { it.name.contains("_${params.golay_index_read}_") },pair.find { it.name.contains('_R2_') })
                    else tuple(sample,pair.find { it.name.contains('_R1_') },pair.find { it.name.contains('_R2_') })
                }
            }
            if(params.demux_mode=='original') {
                golayInput = rows.join(assigned,by:0,failOnDuplicate:true,failOnMismatch:true)
                    .map { sample,json,row,r1,index_read,r2 -> tuple(sample,json,r1,index_read,r2) }
                BUSHMAN_GOLAY_SAMPLE(golayInput,sources,stageHashes.golay)
                ingress = BUSHMAN_GOLAY_SAMPLE.out.reads
            } else ingress = assigned
        } else {
            CREATE_demux_samplesheet_local(NORMALIZE_index_length_local.out.normalized)
            demuxInput = CREATE_demux_samplesheet_local.out.demux_sheet.map { sheet ->
                tuple(file(params.FASTQfolderDir,type:'dir',checkIfExists:true),sheet,params.readStructure)
            }
            DEMUXING_FASTQ_local(demuxInput)
            ingress = DEMUXING_FASTQ_local.out.fastq.flatten()
                .filter { !it.name.startsWith('unmatched.') }
                .map { read -> tuple(read.name.replaceFirst(/\.R[12]\.fq\.gz$/, ''),read) }
                .groupTuple(by:0).map { sample, pair ->
                    if(pair.size()!=2) error "Expected one R1/R2 pair for ${sample}"
                    tuple(sample,pair.find { it.name.endsWith('.R1.fq.gz') },pair.find { it.name.endsWith('.R2.fq.gz') })
                }
        }
    }
    trimInput = rows.join(ingress,by:0,failOnDuplicate:true,failOnMismatch:true)
        .map { sample, json, row, r1, r2 ->
            tuple(sample,json,r1,r2,file(row.vector_fasta,checkIfExists:true),file(row.reference_twobit,checkIfExists:true),row.refGenome)
        }
    BUSHMAN_TRIM(trimInput,sources,blat,stageHashes.trim)
    prepared = BUSHMAN_TRIM.out.prepared.map { sample,row,folder,twobit,genome ->
        def manifest = new groovy.json.JsonSlurper().parseText(folder.resolve('chunks_manifest.json').text)
        def expected = manifest.expected_psl_names as List
        if(manifest.sample != [sample] || expected.toSet().size()!=expected.size())
            error "Invalid chunk manifest for ${sample}"
        tuple(sample,row,folder,twobit,genome,expected)
    }
    completionKeys = prepared.map { sample,row,folder,twobit,genome,expected ->
        tuple(sample,groupKey(sample,expected.size()+1))
    }
    chunks = BUSHMAN_TRIM.out.chunks.join(completionKeys,by:0,failOnDuplicate:true,failOnMismatch:false)
        .flatMap { sample, queries, twobit, genome, key ->
            (queries instanceof List ? queries : [queries]).collect { query -> tuple(sample,query,twobit,genome,key) }
    }
    BUSHMAN_BLAT(chunks,blat)
    // The sentinel completes zero-query samples too. groupKey emits a sample
    // once its expected chunks finish, without waiting for other samples.
    completed = BUSHMAN_BLAT.out.psl.mix(completionKeys.map { sample,key -> tuple(key,null) })
        .groupTuple(by:0).map { key, psls ->
            tuple(key.getGroupTarget(),psls.findAll { it!=null }.sort { it.name })
        }
    callInput = prepared.join(completed,by:0,failOnDuplicate:true,failOnMismatch:true)
        .map { sample,row,folder,twobit,genome,expected,psls ->
            def observed = psls.collect { it.name }
            if(observed.size()!=expected.size() || observed.toSet().size()!=observed.size() || observed.sort()!=expected.sort())
                error "Incomplete or duplicate BLAT chunk set for ${sample}"
            tuple(sample,row,folder,twobit,genome,psls)
        }
    BUSHMAN_CALL(callInput,sources,stageHashes.call)
    pooled = BUSHMAN_CALL.out.calls.map { sample, calls, genome, json ->
        def row = new groovy.json.JsonSlurper().parseText(json.text)
        tuple(genome+'__'+row.Report_Group,calls)
    }.groupTuple(by:0)
    BUSHMAN_POOL(pooled,metadata,sources,stageHashes.pool)
    if(params.generate_plots) BUSHMAN_PLOTS(BUSHMAN_POOL.out.results,metadata,stageHashes.plots)
    emit:
    results = BUSHMAN_POOL.out.results
    provenance = BUSHMAN_VALIDATE.out.provenance
    implementation = BUSHMAN_SOURCE_PROVENANCE.out.provenance
    trim_audit = BUSHMAN_TRIM.out.audit
    calling_audit = BUSHMAN_CALL.out.audit
    pool_audit = BUSHMAN_POOL.out.audit
}
