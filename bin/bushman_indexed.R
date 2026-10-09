# Performance overlay for the frozen, modernized Bushman caller.
# Candidate generation, filters, coordinates, graph and abundance expressions
# remain in vendor/bushman. Replacements below only change lookup/expansion.
bushman_index <- function(values) {
    index <- data.table::data.table(value=as.character(values))
    data.table::setindexv(index, "value")
    index
}

bushman_rows <- function(index, queries) {
    if(!length(queries) || !nrow(index)) return(integer())
    wanted <- data.table::data.table(value=unique(as.character(queries)))
    rows <- index[wanted, on="value", which=TRUE, nomatch=0L]
    sort.int(unique(rows), method="radix")
}

bushman_locus_membership <- function(key_values, queries, locus_ids) {
    if(!length(locus_ids)) return(IRanges::IntegerList())
    index <- bushman_index(key_values)
    loci <- unique(locus_ids)
    first <- match(loci, locus_ids)
    members <- lapply(first, function(i) bushman_rows(index, queries[[i]]))
    # which(%in%) returns ascending original key-row IDs, not query order.
    IRanges::IntegerList(members[match(locus_ids, loci)])
}

bushman_chimera_alignments <- function(hits.R1, hits.R2, reads) {
    if(!nrow(reads)) return(GenomicRanges::GRangesList())
    r1 <- unique(reads$R1); r2 <- unique(reads$R2)
    i1 <- bushman_index(hits.R1$qName); i2 <- bushman_index(hits.R2$qName)
    rows1 <- lapply(r1, function(q) bushman_rows(i1, q))
    rows2 <- lapply(r2, function(q) bushman_rows(i2, q))
    m1 <- match(reads$R1, r1); m2 <- match(reads$R2, r2)
    offset <- length(hits.R2)
    # Keep the original R2-then-R1 order for each original read. Construct the
    # compressed range list in one pass instead of an S4 subset per read.
    rows <- IRanges::IntegerList(lapply(seq_len(nrow(reads)), function(i)
        c(rows2[[m2[i]]], offset + rows1[[m1[i]]])))
    flat <- c(hits.R2, hits.R1)[unlist(rows, use.names=FALSE)]
    names(flat) <- rep(as.character(reads$names), S4Vectors::elementNROWS(rows))
    relist(flat, rows)
}

bushman_multihit_lengths <- function(keys, component_keys) {
    pair_index <- bushman_index(keys$readPairKey)
    id_index <- bushman_index(keys$ID)
    lapply(component_keys, function(component) {
        read_ids <- unique(keys$ID[bushman_rows(pair_index, component)])
        # Preserve the upstream second lookup across ALL rows sharing an ID.
        rows <- bushman_rows(id_index, read_ids)
        data.frame(table(keys$medians[rows]))
    })
}

bushman_group_rows <- function(values) {
    if(!length(values)) return(IRanges::IntegerList())
    runs <- rle(as.character(values))
    if(anyNA(runs$values) || anyDuplicated(runs$values))
        stop("Multihit expansion requires the original sorted nonmissing pair keys")
    ends <- cumsum(runs$lengths)
    starts <- c(1L, head(ends, -1L)+1L)
    result <- IRanges::IntegerList(lapply(seq_along(ends), function(i)
        seq.int(starts[i], ends[i])))
    # The caller retains its original names assignment immediately afterwards.
    result
}

bushman_readpsl_stream <- function(pslFile, toNull=NULL) {
    cols <- c("matches","misMatches","repMatches","nCount","qNumInsert","qBaseInsert",
        "tNumInsert","tBaseInsert","strand","qName","qSize","qStart","qEnd",
        "tName","tSize","tStart","tEnd","blockCount","blockSizes","qStarts","tStarts")
    tables <- lapply(pslFile, function(f) {
        # fread receives the same decompressed text without readLines + paste
        # copies. Retain per-file type inference and input-file ordering.
        tab <- suppressWarnings(data.table::fread(
            cmd=paste("gzip -cd --", shQuote(normalizePath(f, mustWork=TRUE))),
            sep="\t", header=FALSE, showProgress=FALSE, nThread=1L))
        if(!nrow(tab)) return(NULL)
        if(ncol(tab)!=length(cols)) stop("PSL must contain 21 columns: ", f)
        data.table::setnames(tab, cols)
        tab
    })
    result <- data.table::rbindlist(tables)
    if(!nrow(result)) {
        result <- as.data.frame(setNames(lapply(cols, function(n)
            if(n %in% c("strand","qName","tName","blockSizes","qStarts","tStarts"))
                character() else numeric()), cols))
    } else result <- as.data.frame(result)
    if(length(toNull)) result[toNull] <- NULL
    result
}

bushman_optimize <- function(e) {
    data.table::setDTthreads(1L)
    replacements <- list(
        "loci.key$R1.readPairs"=quote(bushman_locus_membership(
            unique_key_pairs$R1, loci.key$R1.qNames, loci.key$R1.loci)),
        "loci.key$R2.readPairs"=quote(bushman_locus_membership(
            unique_key_pairs$R2, loci.key$R2.qNames, loci.key$R2.loci)),
        "chimera.alignments"=quote(bushman_chimera_alignments(hits.R1, hits.R2, chimera.reads)),
        "clusteredMultihitLengths"=quote(bushman_multihit_lengths(multihit.keys, clusteredMultihitNames)),
        "multihit.readPair.read.exp"=quote(bushman_group_rows(multihit.keys$readPairKey)))
    counts <- setNames(integer(length(replacements)), names(replacements))
    rewrite <- function(expr) {
        if(!is.call(expr)) return(expr)
        if(identical(expr[[1]], as.name("<-"))) {
            lhs <- paste(deparse(expr[[2]]), collapse="")
            rhs_head <- if(is.call(expr[[3]])) as.character(expr[[3]][[1]]) else ""
            eligible <- lhs %in% names(replacements)
            if(lhs=="clusteredMultihitLengths") eligible <- eligible && identical(rhs_head,"lapply")
            if(lhs=="multihit.readPair.read.exp") eligible <- eligible && identical(rhs_head,"IntegerList")
            if(eligible) {
                counts[lhs] <<- counts[lhs]+1L
                expr[[3]] <- replacements[[lhs]]
                return(expr)
            }
        }
        # A NULL language argument must remain an argument. [[<- NULL would
        # delete it and corrupt upstream metadata-removal assignments.
        for(i in seq_along(expr)) if(i>1L) expr[i] <- list(rewrite(expr[[i]]))
        expr
    }
    optimized_body <- rewrite(body(e$processAlignments))
    if(!all(counts==1L)) stop("Frozen caller optimization anchors changed: ", paste(counts,collapse=","))
    helper_names <- c("bushman_locus_membership","bushman_chimera_alignments",
        "bushman_multihit_lengths","bushman_group_rows")
    e$kernel_seconds <- list()
    for(name in helper_names) {
        e[[name]] <- local({
            label <- name; kernel <- get(name, mode="function")
            function(...) {
                started <- proc.time()[3]
                result <- kernel(...)
                previous <- e$kernel_seconds[[label]]
                e$kernel_seconds[[label]] <- (if(is.null(previous)) 0 else previous)+unname(proc.time()[3]-started)
                result
            }
        })
    }
    body(e$processAlignments) <- optimized_body
    e$readpsl <- bushman_readpsl_stream
    e$optimization_contract <- list(implementation="indexed", anchors=as.list(counts),
        candidate_generation="unchanged frozen GenomicRanges/BLAT", graph="unchanged frozen igraph",
        abundance="unchanged unique-only Sonic and separate multihits")
    e
}
