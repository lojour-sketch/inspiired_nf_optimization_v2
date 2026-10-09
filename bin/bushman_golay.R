#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly=TRUE)
own <- dirname(normalizePath(sub("^--file=","",grep("^--file=",commandArgs(),value=TRUE)[1])))
source(file.path(own,"bushman_runtime.R"))
e <- bushman_runtime(args[5])
started <- proc.time()
metadata <- if(grepl("\\.json$",args[1])) data.table::as.data.table(jsonlite::fromJSON(args[1])) else
    data.table::fread(args[1],colClasses="character")
idx <- ShortRead::readFastq(args[3])
n_index <- length(idx)
if(n_index) {
    idx <- ShortRead::trimTailw(idx,2,"0",12)
    idx <- idx[width(idx)==max(width(idx))]
    n_quality_retained <- length(idx)
    parts <- split(idx,ceiling(seq_along(idx)/500000))
    for(i in seq_along(parts)) ShortRead::writeFasta(parts[[i]],file=paste0("trimmedI1-",i,".fasta"))
    status <- system2("python3",c(shQuote(file.path(own,"bushman_golay.py")),
        shQuote(file.path(args[5],"intSiteCaller/errorCorrectIndices/golay.py")),
        shQuote(list.files(pattern="^trimmedI1-.*\\.fasta$"))))
    if(status!=0L) stop("Golay correction failed")
    index <- ShortRead::readFasta(list.files(pattern="^correctedI1-.*\\.fasta$"))
} else {
    n_quality_retained <- 0L
    index <- ShortRead::ShortRead(Biostrings::DNAStringSet(),Biostrings::BStringSet())
}
n_decoded <- length(index)
barcodes <- if("bcSeq" %in% names(metadata)) metadata$bcSeq else metadata$index2
index <- index[as.character(sread(index)) %in% barcodes]
ids <- sub(" .*$","",as.character(ShortRead::id(index)))
samples <- metadata$Sample_ID[match(as.character(sread(index)),barcodes)]
counts <- setNames(vapply(barcodes,function(b) sum(as.character(sread(index))==b),integer(1)),metadata$Sample_ID)
for(mate in c("R1","R2")) {
    reads <- ShortRead::readFastq(if(mate=="R1") args[2] else args[4])
    for(i in seq_len(nrow(metadata))) {
        keep <- ids[as.character(sread(index))==barcodes[i]]
        readids <- sub(" .*$","",as.character(ShortRead::id(reads)))
        selected <- match(keep,readids)
        if(anyNA(selected)) stop("Index/read synchronization failure")
        ShortRead::writeFastq(reads[selected],paste0(metadata$Sample_ID[i],".",mate,".fq.gz"),compress=TRUE)
    }
    rm(reads);invisible(gc())
}
jsonlite::write_json(list(input_index_records=n_index,index_quality_full_length_records=n_quality_retained,
    decoded_index_records=n_decoded,assigned_pairs=as.list(counts),
    elapsed_seconds=unname((proc.time()-started)[3]),
    index_rule='Original trimTailw(index, 2, "0", 12), keep maximum length, frozen Golay decode, match bcSeq',
    scope=if(grepl("\\.json$",args[1])) 'One preassigned library; preserve its expected barcode, as historical RUN809' else 'Raw index pool; assign original barcode membership'),
    "golay_metrics.json",pretty=TRUE,auto_unbox=TRUE)
write_runtime("golay_runtime.json")
