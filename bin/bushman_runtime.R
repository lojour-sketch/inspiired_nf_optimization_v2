# Runtime adapters for frozen INSPIIRED sources. Scientific expressions stay
# upstream; changes below repair APIs, input identifiers and empty collections.
bushman_runtime <- function(source_dir, window=5L, reference=NULL, reference_seqinfo=NULL) {
    own <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(), value=TRUE)[1])))
    lib <- file.path(own, "..", "vendor", "bushman_r_library")
    if(dir.exists("/opt/bushman_r_library")) .libPaths(c("/opt/bushman_r_library",.libPaths()))
    else if(dir.exists(lib)) .libPaths(c(normalizePath(lib), .libPaths()))
    packages <- c("plyr","BiocParallel","Biostrings","GenomicAlignments",
        "hiAnnotator","sonicLength","GenomicRanges","ShortRead","igraph","data.table")
    for(p in packages) suppressPackageStartupMessages(library(p,character.only=TRUE))
    # pairwiseAlignment moved to pwalign; scoring arguments remain upstream.
    if(requireNamespace("pwalign",quietly=TRUE)) suppressPackageStartupMessages(library(pwalign))
    e <- new.env(parent=globalenv())
    e$distinct <- dplyr::distinct
    e$count <- plyr::count
    e$arrange <- plyr::arrange
    caller <- file.path(source_dir,"intSiteCaller")
    eval_functions <- function(path) {
        for(expr in parse(path)) {
            if(is.call(expr) && as.character(expr[[1]]) %in% c("<-","=") &&
               is.call(expr[[3]]) && identical(expr[[3]][[1]],as.name("function")))
                eval(expr,e)
        }
    }
    for(name in c("hiReadsProcessor.R","standardization_based_on_clustering.R",
                  "read_psl_files.R","quality_filter.R","intSiteLogic.R"))
        eval_functions(file.path(caller,name))
    # The scheduler/configuration belongs to Nextflow, not the worker function.
    e$source <- function(file,...) {
        if(basename(file)=="programFlow.R") return(invisible(NULL))
        base::source(file,...)
    }
    if(!is.null(reference)) {
        dictionary <- GenomeInfoDb::seqinfo(rtracklayer::TwoBitFile(reference))
        if(!is.null(reference_seqinfo) && nzchar(reference_seqinfo)) {
            frozen_dictionary <- readRDS(reference_seqinfo)
            if(!methods::is(frozen_dictionary,"Seqinfo") ||
               !identical(seqnames(dictionary),seqnames(frozen_dictionary)) ||
               !identical(seqlengths(dictionary),seqlengths(frozen_dictionary)))
                stop("Frozen Seqinfo must exactly match the selected twoBit sequence names, order and lengths")
            dictionary <- frozen_dictionary
        }
        e$get_reference_genome <- function(refGenome) GenomicRanges::GRanges(seqinfo=dictionary)
    }
    # Modern Illumina names need no legacy "read-" prefix. Keep the entire ID.
    text <- paste(deparse(body(e$getTrimmedSeqs),width.cutoff=500L),collapse="\n")
    corrected <- sub('strsplit(z, "-")[[1]][2]', 'z',text,fixed=TRUE)
    if(identical(corrected,text)) stop("Modern-header compatibility adapter did not match the frozen function")
    text <- corrected
    for(mate in c("reads.p","reads.l")) text <- sub(
        paste0("base::subset(",mate,", width(",mate,") > mingDNA)"),
        paste0(mate,"[width(",mate,") > mingDNA]"),text,fixed=TRUE)
    # Empty inputs must not leave uninitialized sequence objects.
    text <- sub('seqs <- x\\[\\[1\\]\\]', paste0('seqs <- x[[1]]\n',
        'trimmedSeqs <- DNAStringSet(); trimmedqSeqs <- BStringSet()'),text)
    body(e$getTrimmedSeqs) <- parse(text=text)[[1]]
    # Retain the original two-flank algorithm; add the known leading-UMI
    # library layout adapter only when the first linker symbol is N.
    original_linker <- e$trim_primerIDlinker_side_reads
    e$trim_primerIDlinker_side_reads <- function(reads.l,linker,maxMisMatch=3) {
        pos <- unlist(gregexpr("N",linker))
        if(min(pos)!=1L) return(original_linker(reads.l,linker,maxMisMatch))
        count <- max(pos)-min(pos)+1L; link2 <- substr(linker,max(pos)+1L,nchar(linker))
        aln <- pairwiseAlignment(pattern=subseq(reads.l,max(1,count-1L),count+nchar(link2)+1L),
            subject=link2,substitutionMatrix=nucleotideSubstitutionMatrix(match=1,mismatch=0,baseOnly=TRUE),
            gapOpening=0,gapExtension=1,type="overlap")
        df <- e$PairwiseAlignmentsSingleSubject2DF(aln,max(1,count-1L)-1L)
        keep <- df$score>=nchar(link2)-maxMisMatch
        list(reads.l=subseq(reads.l[keep],df$end[keep]+1L),
             primerID=subseq(reads.l[keep],1,count))
    }
    # GRanges is already a range vector: modern unlist(GRanges) is invalid.
    txt <- paste(deparse(body(e$dereplicateSites),width.cutoff=500L),collapse="\n")
    txt <- sub('unlist\\(reduce\\(sites.reduced, min.gapwidth = 0L, with.revmap = TRUE\\)\\)',
               'reduce(sites.reduced, min.gapwidth = 0L, with.revmap = TRUE)',txt)
    body(e$dereplicateSites) <- parse(text=txt)[[1]]
    # Preserve the original <=window assertion, frequency, tie and one-pass
    # rules. Explicit self overlaps replace the obsolete one-argument API.
    e$self_overlaps <- function(query,ignoreSelf=TRUE,ignoreRedundant=FALSE,select="all",maxgap) {
        hits <- GenomicRanges::findOverlaps(query,query,maxgap=maxgap-1L,select=select)
        hits[S4Vectors::queryHits(hits)!=S4Vectors::subjectHits(hits)]
    }
    txt <- paste(deparse(body(e$clusterSites),width.cutoff=500L),collapse="\n")
    txt <- sub('findOverlaps\\(sites.gr, ignoreSelf', 'self_overlaps(sites.gr, ignoreSelf',txt)
    body(e$clusterSites) <- parse(text=txt)[[1]]
    formals(e$clusterSites)$windowSize <- as.integer(window)
    # Iterating seq_len(0) fixes legacy 1:0 loops, without altering nonempty data.
    rewrite <- function(x) {
        if(!is.call(x)) return(x)
        if(identical(x[[1]],as.name(":")) && identical(x[[2]],1L) &&
           is.call(x[[3]]) && as.character(x[[3]][[1]]) %in% c("length","nrow"))
            return(as.call(list(as.name("seq_len"),rewrite(x[[3]]))))
        if(identical(x[[1]],as.name(":")) && identical(x[[2]],1) &&
           is.call(x[[3]]) && as.character(x[[3]][[1]]) %in% c("length","nrow"))
            return(as.call(list(as.name("seq_len"),rewrite(x[[3]]))))
        if(identical(x[[1]],as.name("::")) && identical(x[[2]],as.name("igraph"))) {
            if(identical(x[[3]],as.name("clusters"))) x[[3]] <- as.name("components")
            if(identical(x[[3]],as.name("graph.edgelist"))) x[[3]] <- as.name("graph_from_edgelist")
        }
        for(i in seq_along(x)) if(i>1L) x[[i]] <- rewrite(x[[i]])
        x
    }
    body(e$processAlignments) <- rewrite(body(e$processAlignments))
    txt <- paste(deparse(body(e$processAlignments),width.cutoff=500L),collapse="\n")
    txt <- sub("seq_len(length(chimera.reads))","seq_len(nrow(chimera.reads))",txt,fixed=TRUE)
    body(e$processAlignments) <- parse(text=txt)[[1]]
    # Record stage timings around original statements; no filtering shortcut.
    e$timings <- list(); e$last_phase <- NULL; e$last_time <- proc.time()[3]
    e$mark <- function(label) {
        now <- proc.time()[3]
        if(!is.null(e$last_phase)) e$timings[[e$last_phase]] <- (if(is.null(e$timings[[e$last_phase]])) 0 else e$timings[[e$last_phase]]) + unname(now-e$last_time)
        e$last_phase <- label;e$last_time <- now
        invisible(NULL)
    }
    instrument <- function(fun,labels) {
        statements <- as.list(body(fun));out <- statements[1]
        for(statement in statements[-1]) {
            if(is.call(statement) && identical(statement[[1]],as.name("<-"))) {
                name <- paste(deparse(statement[[2]]),collapse="")
                if(name %in% names(labels)) out <- c(out,list(substitute(mark(LABEL),list(LABEL=labels[[name]]))))
            }
            out <- c(out,list(statement))
        }
        body(fun) <- as.call(out);fun
    }
    e$getTrimmedSeqs <- instrument(e$getTrimmedSeqs,c(reads="read_loading",
        reads.p="LTR_or_readthrough",readslprimer="linker_processing",
        vqName="vector_BLAT",reads.p.u="sequence_dereplication"))
    e$processAlignments <- instrument(e$processAlignments,c(psl.R2="PSL_loading",
        unique_key_pairs="locus_pairing",failedReads="chimera_expansion",
        uniq.read.loci.mat="unique_read_expansion",unclusteredMultihits="multihit_graph_and_expansion"))
    # Typed empty PSL tables and explicit cmd= replace fread's shell inference.
    e$readpsl <- function(pslFile,toNull=NULL) {
        cols <- c("matches","misMatches","repMatches","nCount","qNumInsert","qBaseInsert",
            "tNumInsert","tBaseInsert","strand","qName","qSize","qStart","qEnd",
            "tName","tSize","tStart","tEnd","blockCount","blockSizes","qStarts","tStarts")
        tables <- lapply(pslFile,function(f) {
            lines <- readLines(gzfile(f),warn=FALSE)
            if(!length(lines)) return(NULL)
            data.table::fread(text=paste(lines,collapse="\n"),sep="\t",header=FALSE,col.names=cols)
        })
        result <- data.table::rbindlist(tables)
        if(!nrow(result)) {
            result <- as.data.frame(setNames(lapply(cols,function(n)
                if(n %in% c("strand","qName","tName","blockSizes","qStarts","tStarts")) character() else numeric()),cols))
        } else result <- as.data.frame(result)
        if(length(toNull)) result[toNull] <- NULL
        result
    }
    e
}

write_runtime <- function(path) {
    versions <- vapply(loadedNamespaces(),function(p) as.character(packageVersion(p)),character(1))
    jsonlite::write_json(list(R=R.version.string,packages=as.list(versions),
        implementation_sha256=Sys.getenv("BUSHMAN_IMPLEMENTATION_SHA256",unset=NA_character_),
        adapters=c("namespace_qualification","modern_headers","leading_UMI_layout","empty_collections","XStringSet_subset_API","GRanges_reduce_API",
            "explicit_self_overlaps_and_inclusive_window","igraph_API","chimera_archive_row_index",
            "Nextflow_scheduler_boundary")),
        path,pretty=TRUE,auto_unbox=TRUE)
}
