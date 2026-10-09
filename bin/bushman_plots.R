#!/usr/bin/env Rscript
# Presentation only: consume archived caller/report outputs; never refit,
# merge sites, discard short fragments, or overwrite scientific inputs.
args <- commandArgs(TRUE)
if(!length(args) %in% c(5L,6L)) stop("Usage: REPORT_DIRECTORY METADATA WINDOW_BP OUTPUT_DIRECTORY ANALYSIS_LABEL [VALIDATION_LABEL]")
suppressPackageStartupMessages(library(GenomicRanges))
suppressPackageStartupMessages(library(data.table))
report <- normalizePath(args[1]); meta <- fread(args[2], colClasses="character")
out <- args[4]; dir.create(out,recursive=TRUE,showWarnings=FALSE)
label <- args[5]; validation <- if(length(args)==6L) args[6] else "Validate against the original pipeline before biological interpretation."
load_table <- function(name) fread(file.path(report,name),data.table=FALSE)
abundance <- load_table("original_unique_abundance.tsv")
status <- load_table("abundance_status.tsv")
multihits <- load_table("original_multihit_abundance.tsv")
pooled <- readRDS(file.path(report,"pooled_sites.rds"))$standardized
if(length(pooled)) {
    points <- flank(pooled,-1,start=TRUE)
    geometry <- data.frame(GTSP=as.character(points$GTSP),seqnames=as.character(seqnames(points)),
        start=start(points),strand=as.character(strand(points)),posid=as.character(points$posid))
    geometry <- geometry[!duplicated(geometry),,drop=FALSE]
} else geometry <- data.frame(GTSP=character(),seqnames=character(),start=integer(),strand=character(),posid=character())
specimens <- unique(meta$Specimen_ID[meta$Specimen_ID %in% status$GTSP])
if(!length(specimens)) stop("No report specimens match the metadata")
esc <- function(x) {
    x <- gsub("&","&amp;",as.character(x),fixed=TRUE)
    x <- gsub("<","&lt;",x,fixed=TRUE); gsub('"',"&quot;",x,fixed=TRUE)
}
empty_plot <- function(title,text) { plot.new();title(main=title);text(.5,.5,text,cex=.9) }
caption <- function() mtext(paste(label, "|",validation),side=1,line=4,cex=.55)
chromosomes <- function(x) {
    canonical <- c(paste0("chr",1:22),"chrX","chrY","chrM")
    c(intersect(canonical,unique(x)),sort(setdiff(unique(x),canonical)))
}
summaries <- list(); links <- character(); files <- character()
for(specimen in specimens) {
    stopifnot(grepl("^[A-Za-z0-9][A-Za-z0-9_.-]*$",specimen))
    a <- abundance[abundance$GTSP==specimen,,drop=FALSE]
    g <- geometry[geometry$GTSP==specimen,,drop=FALSE]
    st <- status[status$GTSP==specimen,,drop=FALSE]
    valid <- nrow(a)>0L && st$status[1] %in% c("estimated","estimated_with_warnings")
    if(valid) {
        stopifnot(all(is.finite(a$estAbund)),all(a$estAbund>=0),sum(a$estAbund)>0)
        a <- a[order(-a$estAbund,a$posid),,drop=FALSE]
        p <- a$estAbund/sum(a$estAbund)
        stopifnot(isTRUE(all.equal(p,as.numeric(a$estAbundProp),tolerance=1e-10)))
    } else p <- numeric()
    rows <- meta[meta$Specimen_ID==specimen,,drop=FALSE]
    expected <- if("Expected_clonal_fraction" %in% names(rows)) unique(rows$Expected_clonal_fraction) else NA_character_
    if(length(expected)!=1L) stop("Conflicting Expected_clonal_fraction metadata")
    summary <- data.frame(specimen=specimen,abundance_status=st$status[1],
        standardized_unique_sites=nrow(g),abundance_sites=if(valid)nrow(a) else NA_integer_,
        estimated_unique_site_abundance=if(valid)sum(a$estAbund) else NA_real_,
        largest_unique_site_fraction=if(valid)max(p) else NA_real_,
        top10_unique_site_fraction=if(valid)sum(head(p,10)) else NA_real_,
        shannon_unique_site_diversity=if(valid)-sum(p[p>0]*log(p[p>0])) else NA_real_,
        simpson_unique_site_concentration=if(valid)sum(p^2) else NA_real_,
        effective_unique_sites=if(valid)1/sum(p^2) else NA_real_,
        expected_clonal_cell_fraction=expected,site_window_bp=as.integer(args[3]),
        abundance_method=if(valid)a$abundance_method[1] else "unavailable",
        analysis_label=label,validation_label=validation)
    summaries[[specimen]] <- summary
    table_name <- paste0("annotated_points_",specimen,".tsv")
    if(valid) fwrite(a,file.path(out,table_name),sep="\t",na="NA")
    else fwrite(g,file.path(out,table_name),sep="\t",na="NA")
    fwrite(st,file.path(out,paste0("abundance_status_",specimen,".tsv")),sep="\t")
    # Location/annotation figures analogous to legacy Step 17, using the
    # pipeline's frozen Bushman RefSeq annotations rather than another database.
    figure <- paste0("fig_points_",specimen,".pdf")
    pdf(file.path(out,figure),width=11,height=7,onefile=TRUE)
    par(mar=c(7,5,4,2))
    chr <- chromosomes(g$seqnames)
    if(nrow(g)) {
        plot(g$start/1e6,match(g$seqnames,chr),pch=16,col="#2563ab80",cex=.65,
            yaxt="n",xlab="Insertion coordinate (Mb; 1-based)",ylab="Chromosome",
            main=paste("Insertion locations:",specimen),ylim=c(.5,length(chr)+.5))
        axis(2,at=seq_along(chr),labels=chr,las=2,cex.axis=.6);caption()
        tab <- table(factor(g$seqnames,levels=chr))
        barplot(tab,las=2,cex.names=.6,col="#2563ab",ylab="Unique standardized sites",
            main="Chromosome distribution (site counts)");caption()
    } else empty_plot("Insertion locations","No unique sites")
    if(valid && "inGene" %in% names(a)) {
        category <- ifelse(is.na(a$inGene)|a$inGene=="","Unavailable",ifelse(a$inGene=="FALSE","Intergenic","Intragenic"))
        tab <- table(category)
        pie(tab,labels=paste(names(tab),tab,sep="\n"),col=c("#2563ab","#f59e0b","#aaa")[seq_along(tab)],
            main="Gene overlap: frozen Bushman RefSeq (unweighted sites)");caption()
        if("X5pnearest_refSeq_geneDist" %in% names(a)) {
            distance <- a$X5pnearest_refSeq_geneDist
            if(any(is.finite(distance))) {
                hist(sign(distance)*log10(1+abs(distance)),breaks=40,col="#2563ab",border="white",
                    xlab="sign(distance) × log10(1 + |distance in bp|)",
                    main="Distance to nearest 5-prime RefSeq feature (upstream annotation)");caption()
            }
        }
    } else empty_plot("Gene annotation",if(valid)"Annotation disabled" else "Sonic abundance unavailable; annotation report unavailable")
    dev.off()
    clonal <- paste0("clonality_",specimen,".pdf")
    pdf(file.path(out,clonal),width=11,height=7,onefile=TRUE);par(mar=c(7,5,4,2))
    if(valid) {
        plot(seq_along(p),100*p,log="x",pch=16,cex=.5,col="#2563ab",
            xlab="Unique insertion-site rank (log scale)",ylab="Unique-site Sonic abundance (%)",
            main=paste("Abundance rank:",specimen));caption()
        plot(c(0,seq_along(p)),100*c(0,cumsum(p)),type="l",lwd=2,col="#2563ab",
            xlab="Highest-abundance unique sites included",ylab="Cumulative unique-site abundance (%)",
            main="Cumulative abundance: dominant sites first",ylim=c(0,100));caption()
        top <- head(a,20); heights <- 100*head(p,20)
        barplot(heights,names.arg=top$posid,las=2,cex.names=.6,col="#2563ab",
            ylab="Unique-site Sonic abundance (%)",main="Top 20 insertion sites");caption()
        remaining <- max(0,1-sum(head(p,10)))
        pie(c(head(p,10),remaining),labels=c(paste0("Rank ",seq_len(min(10,length(p)))),"Other unique sites"),
            main="Unique-site abundance shares (multihits excluded)");caption()
    } else empty_plot(paste("Clonality:",specimen),paste("Abundance unavailable\n",st$warning[1],
        "\nOriginal estimator inputs retained. No substitute estimate."))
    m <- multihits[multihits$GTSP==specimen,,drop=FALSE]
    if(nrow(m)) {
        cluster <- m[!duplicated(m$multihitID),,drop=FALSE]
        cluster <- cluster[order(-cluster$estAbund),,drop=FALSE]
        barplot(head(cluster$estAbund,20),names.arg=head(cluster$multihitID,20),las=2,cex.names=.6,
            col="#d97706",ylab="Distinct fragment lengths (original multihit metric)",
            main="Multihit clusters: separate from unique-site Sonic abundance");caption()
    }
    dev.off()
    links <- c(links,paste0('<li><strong>',esc(specimen),'</strong>: ',esc(st$status[1]),
        ' — <a href="',figure,'">location and annotation PDF</a>; <a href="',clonal,
        '">clonality PDF</a>; <a href="',table_name,'">site table</a></li>'))
    files <- c(files,figure,clonal,table_name)
}
summary_table <- rbindlist(summaries);fwrite(summary_table,file.path(out,"clonality_summary.tsv"),sep="\t",na="NA")
html <- paste0('<!doctype html><html><head><meta charset="utf-8"><title>Insertion-site and clonality results</title>',
    '<style>body{font:17px system-ui;max-width:1050px;margin:40px auto;padding:0 24px;line-height:1.55}li{margin:12px 0}.notice{padding:16px;background:#fff3cd}a{color:#165a9e}</style></head><body>',
    '<h1>Insertion-site and clonality results</h1><p>',esc(label),'</p><p class="notice">',esc(validation),'</p>',
    '<p>Coordinate window: ',esc(args[3]),' bp. Coordinates are 1-based, strand-aware insertion points. ',
    'Figures use the scientific outputs already generated; plotting does not change them.</p>',
    '<p>Rank, cumulative and top-site plots use rounded unique-site Sonic abundance normalized within the specimen. ',
    'Multihit clusters retain the original distinct-fragment-length metric and are shown separately. ',
    'Unique-site fractions and diversity describe insertion sites; they do not identify cellular clones or directly measure clonal cell fractions. ',
    'Gene overlap and nearest 5-prime feature plots use the frozen Bushman RefSeq annotation. They are not ChIPseeker exon/promoter classifications.</p>',
    '<p>Expected mixture fractions, when supplied as Expected_clonal_fraction metadata, refer to cells and are reported for context. ',
    'Unavailable abundance is NA, not zero; no invalid short fragment is filtered for plotting.</p>',
    '<p><a href="clonality_summary.tsv">All-specimen metrics</a></p><ul>',paste(links,collapse="\n"),'</ul></body></html>')
writeLines(html,file.path(out,"index.html"))
jsonlite::write_json(list(state="completed",analysis_label=label,validation_label=validation,
    site_window_bp=as.integer(args[3]),specimens=specimens,figures=files,
    inputs_modified=FALSE,abundance_refit=FALSE,unavailable_status_retained=TRUE),
    file.path(out,"plot_manifest.json"),pretty=TRUE,auto_unbox=TRUE)
cat("PASS: location, annotation and clonality reports for",length(specimens),"specimens\n")
