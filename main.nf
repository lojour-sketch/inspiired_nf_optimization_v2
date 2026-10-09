nextflow.enable.dsl=2
params.BCLorFASTQ = ''
params.samplesheet = ''
params.runfolderDir = ''
params.FASTQfolderDir = ''
params.fastq_manifest = ''
params.readStructure = ''
params.instrument = ''
params.projectName = ''
params.outdir = ''
params.reference_manifest = ''
params.demux_mode = 'original'
params.annotation_mode = 'bushman'
params.abundance_method = 'sonic'
params.abundance_failure_policy = 'fail'
params.site_window_bp = 5
params.blat_chunk_size = 30000
params.caller_impl = 'compact'
params.generate_plots = true
params.blat_binary = "${projectDir}/vendor/bushman/tools/blat"
params.bcl_bases_mask = 'I20Y159,I12,Y143'
params.golay_bcl_bases_mask = 'Y179,I12,Y143'
params.golay_index_read = 'I2'
params.bushman_r_container = "${projectDir}/containers/bushman_runtime.sif"
if(!params.samplesheet || !params.projectName || !params.outdir || !params.runfolderDir)
    error 'Provide --samplesheet, --runfolderDir, --projectName and --outdir'
if(!(params.BCLorFASTQ in ['BCL','FASTQ'])) error '--BCLorFASTQ must be BCL or FASTQ'
if(!params.reference_manifest) error 'Provide --reference_manifest selecting the exact original-compatible twoBit and frozen RefSeq annotation'
if(!(params.demux_mode in ['original','fixed','golay'])) error '--demux_mode must be original, fixed or golay'
if(params.demux_mode=='golay' && params.BCLorFASTQ=='BCL') error 'The original Golay adapter accepts raw FASTQs; use the fixed BCL input adapter for inline dual-index BCL libraries'
if(!(params.golay_index_read in ['I1','I2'])) error '--golay_index_read must be I1 or I2'
if(!(params.annotation_mode in ['bushman','off'])) error '--annotation_mode must be bushman or off'
if(!(params.abundance_method in ['sonic','fragment_diversity'])) error '--abundance_method must be sonic or fragment_diversity'
if(!(params.abundance_failure_policy in ['fail','record'])) error '--abundance_failure_policy must be fail or record'
if(!(params.site_window_bp.toString()==~/[0-9]+/)) error '--site_window_bp must be a nonnegative integer'
if(!(params.blat_chunk_size.toString()==~/[1-9][0-9]*/)) error '--blat_chunk_size must be a positive integer'
if(!(params.caller_impl in ['compact','indexed','original'])) error '--caller_impl must be compact, indexed or original'
if(params.BCLorFASTQ=='FASTQ' && !params.fastq_manifest) {
    if(!params.FASTQfolderDir) error 'Provide --FASTQfolderDir or --fastq_manifest'
    if(params.demux_mode=='fixed' && (!params.readStructure || !(params.instrument in ['MiSeq','NextSeq2000','NextSeq500'])))
        error 'The fixed FASTQ input adapter requires --readStructure and --instrument'
}
include { ORIGINAL_COMPATIBLE } from './subworkflows/local/bushman/main'
workflow { ORIGINAL_COMPATIBLE() }
