# Experimental execution overlay: frozen candidate generation and science stay
# in intSiteLogic.R. Store membership once per locus, not once per candidate.
bushman_load_pairing <- function(source_path, build_root=file.path(getwd(),".bushman_native")) {
    source_path <- normalizePath(source_path, mustWork=TRUE)
    fingerprint <- unname(tools::md5sum(source_path))
    directory <- file.path(build_root,paste0(fingerprint,"-R",getRversion()))
    dir.create(directory,recursive=TRUE,showWarnings=FALSE)
    directory <- normalizePath(directory,mustWork=TRUE)
    library <- file.path(directory,paste0("bushman_pairing",.Platform$dynlib.ext))
    if(!file.exists(library)) {
        file.copy(source_path,file.path(directory,"bushman_pairing.c"),overwrite=TRUE)
        old <- getwd(); on.exit(setwd(old),add=TRUE); setwd(directory)
        status <- system2(file.path(R.home("bin"),"R"),
            c("CMD","SHLIB","bushman_pairing.c"),stdout="compile.log",stderr="compile.log")
        if(status!=0L || !file.exists(library))
            stop("Compact caller compilation failed; inspect ",file.path(directory,"compile.log"))
        setwd(old)
    }
    dll <- dyn.load(library)
    list(symbol=getNativeSymbolInfo("bushman_intersections",dll),
        provenance=list(source_md5=fingerprint,library_md5=unname(tools::md5sum(library)),
            R_version=R.version.string,platform=R.version$platform,
            build_directory=directory,build_command="R CMD SHLIB bushman_pairing.c"))
}

bushman_compact_members <- function(key_values, hits, reduced, locus_ids, block_size=4096L) {
    loci <- unique(as.integer(locus_ids))
    if(!length(loci)) return(list(values=integer(),ends=integer(),slots=integer()))
    index <- data.table::data.table(value=as.character(key_values),rowid=seq_along(key_values))
    data.table::setindexv(index,"value")
    blocks <- seq.int(1L,length(loci),by=block_size)
    values <- lengths <- vector("list",length(blocks))
    for(block in seq_along(blocks)) {
        positions <- seq.int(blocks[block],min(length(loci),blocks[block]+block_size-1L))
        revmap <- reduced$revmap[loci[positions]]
        queries <- data.table::data.table(
            # Keep the original integer coercion before %in% matching, including NA.
            value=as.character(as.integer(hits$qName[unlist(revmap,use.names=FALSE)])),
            slot=rep(seq_along(positions),S4Vectors::elementNROWS(revmap)))
        queries <- unique(queries)
        rows <- index[queries,on="value",allow.cartesian=TRUE,nomatch=0L,
            .(slot=i.slot,rowid=x.rowid)]
        rows <- unique(rows)
        data.table::setorderv(rows,c("slot","rowid"))
        values[[block]] <- rows$rowid
        lengths[[block]] <- tabulate(rows$slot,nbins=length(positions))
    }
    ends <- cumsum(unlist(lengths,use.names=FALSE))
    if(length(ends) && tail(ends,1)> .Machine$integer.max)
        stop("Membership exceeds IRanges integer partition limit; no data truncated")
    list(values=as.integer(unlist(values,use.names=FALSE)),ends=as.integer(ends),
        slots=match(locus_ids,loci))
}

bushman_compact_support <- function(keys,hits.R1,hits.R2,red.R1,red.R2,loci.key,native,record) {
    started <- proc.time()[3]
    a <- bushman_compact_members(keys$R1,hits.R1,red.R1,loci.key$R1.loci)
    b <- bushman_compact_members(keys$R2,hits.R2,red.R2,loci.key$R2.loci)
    built <- proc.time()[3]
    support <- .Call(native$symbol,a$values,a$ends,b$values,b$ends,a$slots,b$slots)
    intersected <- proc.time()[3]
    result <- relist(as.character(keys$readPairKey[support$values]),
        IRanges::PartitioningByEnd(support$ends))
    record(list(membership_seconds=unname(built-started),
        intersection_seconds=unname(intersected-built),
        assembly_seconds=unname(proc.time()[3]-intersected),
        candidate_pairs=nrow(loci.key),distinct_R1_loci=length(a$ends),distinct_R2_loci=length(b$ends),
        R1_membership_rows=length(a$values),R2_membership_rows=length(b$values),
        support_rows=length(support$values)))
    result
}

bushman_compact_optimize <- function(e,native) {
    e <- bushman_optimize(e) # retain the independently validated expansion/PSL overlays
    replacements <- list(
        "R1.loci"=quote(GenomicRanges::granges(red.hits.R1,use.mcols=FALSE)[queryHits(pairs)]),
        "R2.loci"=quote(GenomicRanges::granges(red.hits.R2,use.mcols=FALSE)[subjectHits(pairs)]),
        "loci.key$R1.qNames"=quote(NULL),"loci.key$R2.qNames"=quote(NULL),
        "loci.key$R1.readPairs"=quote(NULL),"loci.key$R2.readPairs"=quote(NULL),
        "paired.loci$readPairKeys"=quote(bushman_compact_support(unique_key_pairs,
            hits.R1,hits.R2,red.hits.R1,red.hits.R2,loci.key)))
    counts <- setNames(integer(length(replacements)),names(replacements))
    rewrite <- function(expr) {
        if(!is.call(expr)) return(expr)
        if(identical(expr[[1]],as.name("<-"))) {
            lhs <- paste(deparse(expr[[2]]),collapse="")
            eligible <- lhs %in% names(replacements)
            # Later R1/R2.loci filtering assignments stay in the frozen body.
            if(lhs=="R1.loci") eligible <- eligible && identical(expr[[3]],quote(red.hits.R1[queryHits(pairs)]))
            if(lhs=="R2.loci") eligible <- eligible && identical(expr[[3]],quote(red.hits.R2[subjectHits(pairs)]))
            if(eligible) {
                counts[lhs] <<- counts[lhs]+1L
                expr[3] <- list(replacements[[lhs]])
                return(expr)
            }
        }
        for(i in seq_along(expr)) if(i>1L) expr[i] <- list(rewrite(expr[[i]]))
        expr
    }
    updated <- rewrite(body(e$processAlignments))
    if(!all(counts==1L)) stop("Frozen compact pairing anchors changed: ",paste(counts,collapse=","))
    e$compact_pairing_metrics <- list()
    kernel <- bushman_compact_support
    e$bushman_compact_support <- function(...) kernel(...,native=native,
        record=function(metrics) e$compact_pairing_metrics <- metrics)
    body(e$processAlignments) <- updated
    e$optimization_contract$implementation <- "compact"
    e$optimization_contract$compact_anchors <- as.list(counts)
    e$optimization_contract$pairing <- "one membership per distinct locus; sorted original row-ID intersections; unchanged pair ordering"
    e$optimization_contract$native <- native$provenance
    e
}
