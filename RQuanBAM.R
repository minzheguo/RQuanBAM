

#' BAM Files to Counts and Normalized Counts (CPM and simplified RPKM/TPM)
#' 
#' require R environment 3.2.1
#' inconsistencies may exist if older versions of R is used, 
#' causing the error: could not find symbol "keepNA" in environment of the generic function
#' 
#' @param in.dir (character) the full name of the directory that contains all the input files, including bam files and a gtf/gff file
#' @param gff.format (character) the format of gene annotation file: gtf or gff3
#' @param bam.yieldSize (numeric) Number of records to yield each time the file is read from with scanBam, see ?Rsamtools::BamFileList
#' @param chrominfo (data.frame/seqinfo) Data frame containing information about the chromosomes, see ?GenomicFeatures::makeTxDbFromGFF
#' @param feature.type (character) One of "gene", "exon", "cds" or "tx". Determines the grouping, see ?GenomicFeatures::exonsBy
#' @param mode (character) mode can be one of the pre-defined count methods such as "Union", "IntersectionStrict", or "IntersectionNotEmpty" or it a user supplied count function., see ?GenomicAlignments::summarizeOverlaps
#' @param singleEnd (logical) (Default TRUE) A logical indicating if reads are single or paired-end, see ?GenomicAlignments::summarizeOverlaps
#' @param ignore.strand (logical) A logical indicating if strand should be considered when matching, see ?GenomicAlignments::summarizeOverlaps
#' @param fragments (logical) (Default FALSE) A logical. applied to paired-end data only. see ?GenomicAlignments::summarizeOverlaps
#' @param norm.method (character) the method for counts normalization: "none", "TMM", "upperquartile", "RLE", "DESeq2". see ?edgeR::calcNormFactors    
#' @param percentage (numeric) the percentage of samples that will be quantified
#' @param parallel (logical) A logical indicating whether to use parallel evaluation 
#' @return a list of two elements: counts (rawCounts) and cpm  
#' @export
quantifyFromBam <- function(in.dir,  
                    gff.format="gtf", 
                    bam.yieldSize=2000000, 
                    chrominfo = NULL, 
                    feature.type="gene",
                    mode="union", 
                    singleEnd=TRUE, 
                    ignore.strand=TRUE, 
                    fragments=FALSE,
                    norm.method="none", 
                    percentage=1,
                    parallel=TRUE, export=T, verbose=T) {
    
    # input validation
    gfffs <- c("auto", "gtf", "gff3")
    gff.id <- pmatch(gff.format, gfffs)
    if (is.na(gff.id)) {
        stop("invalid gff format. should be auto, gtf or gff3")
    }
    gff.format <- gfffs[gff.id]
    
    norm.methods <- c("none", "TMM", "upperquartile", "RLE", "DESeq2")
    norm.id <- pmatch(norm.method, norm.methods)
    if (is.na(norm.id)) {
        stop("invalid normalization method. should be none, TMM, upperquantile, or RLE")
    }
    norm.method <- norm.methods[norm.id]
    
    # resolving dependencies
    if (!require("Rsamtools")) {
        source("http://bioconductor.org/biocLite.R")
        biocLite("Rsamtools")
    }
    if (!require("GenomicFeatures")) {
        source("http://bioconductor.org/biocLite.R")
        biocLite("GenomicFeatures")
    }
    if (!require("GenomicFeatures")) {
        source("http://bioconductor.org/biocLite.R")
        biocLite("GenomicFeatures")
    }
    if (!require("GenomicAlignments")) {
        source("http://bioconductor.org/biocLite.R")
        biocLite("GenomicAlignments")
    }
    if (!require("edgeR")) {
        source("http://bioconductor.org/biocLite.R")
        biocLite("edgeR")
    }
    
    # Indicate in Bioconductor that these files are BAM files using the BamFileList function. 
    # Here we also specify details about how the BAM files should be treated, e.g., only process 2 million reads at a time. 
    # See ?BamFileList for more information.
    # require("Rsamtools")
    
    bam.filenames <- list.files(in.dir, pattern="\\.bam$", full.name=TRUE)
    if (length(bam.filenames) < 1) {
        stop("No bam files detected.")
    } else {
        if (verbose) {
            cat("\t", length(bam.filenames), " bam files detected\n", sep="")
        }
    }
    
    gff.filename <- list.files(in.dir, pattern="\\.gtf$|\\.gff3$", full.name=TRUE)
    if (length(gff.filename) > 1) {
        stop("Multiple annotation files detected.")
    } else if (length(gff.filename)==0) {
        stop("No annotation file detected.")
    } else {
        if (verbose) {
            cat("\tAnnotation file detected: ", gff.filename, "\n", sep="")
        }
    }
    
    if (verbose) {
        cat("\tDefining bam file list ... ")
    }
    bamfiles <- BamFileList(bam.filenames, yieldSize=bam.yieldSize)
    if (verbose) {
        cat("Ok\n")
    }
    
    sis <- lapply(bamfiles, seqinfo)
    bam.chrom.names <- lapply(sis, names)
    
    # definining gene models from gfffile
    # require("GenomicFeatures")
    
    if (verbose) {
        cat("\tDefining gene models ...")
    }
    if (is.null(chrominfo)) {
        txdb <- makeTxDbFromGFF(gff.filename, format=gff.format, circ_seqs=character())
    } else {
        txdb <- makeTxDbFromGFF(gff.filename, format=gff.format, circ_seqs=character(), chrominfo=chrominfo)
    }
    if (verbose) {
        cat("Ok\n")
    }
    
    cat("\tChecking the consistencies of chromosome names ...")
    txdb.chrom.names <- names(seqinfo(txdb))
    chrom.names.common <- lapply(bam.chrom.names, intersect, txdb.chrom.names)
    chrom.names.common.len <- lapply(chrom.names.common, length)
    
    if (any(chrom.names.common.len==0)) {
        stop("Inconsistencies in chromosom names exist between bam files and the annotation file.")
    }
    if (verbose) {
        cat("Ok\n")
    }
    
    if (verbose) {
        cat("\tConstructing", feature.type, "annotations ...")
    }
    ebg <- exonsBy(txdb, by=feature.type)
    if (verbose) {
        cat("Ok\n")
    }
    
    # for multi-core computation
    # require("BiocParallel")
    if (parallel) {
        if (!require("BiocParallel")) {
            source("http://bioconductor.org/biocLite.R")
            biocLite("BiocParallel")
        }
        register(SerialParam())
        cat("\tParallel processing enabled\n")
    }
    
    # reads counting
    
    n <- max(floor(length(bamfiles)*percentage), 1)
    if (verbose) {
        cat("\tReads counting for ", n, " out of ", length(bamfiles), " samples", sep="")
    }
    
    
    se <- summarizeOverlaps(features=ebg, reads=bamfiles[1:n], mode="Union", 
                            singleEnd=singleEnd, ignore.strand=ignore.strand, fragments=fragments )
    
    counts.raw <- assay(se)
    rm(se)
    
    
    if (verbose) {
        cat("done\n")
    }
    
    # the library sizes will be computed from the column sums of the counts
    dge <- DGEList(counts=counts.raw)
    
    if (verbose) {
        cat("\tNormalizing counts using ", norm.method, " method...", sep="")
    }
    # normalization
    if (norm.method=="DESeq2") {
        
        if (!require("DESeq2")) {
            source("http://bioconductor.org/biocLite.R")
            biocLite("DESeq2")
        }
        
        sf <- NULL
        sf <- estimateSizeFactorsForMatrix(counts.raw)
        counts.norm <- t(t(counts.raw)/sf)
        dge <- DGEList(counts=counts.norm)
    } else if (norm.method=="none") {
    } else {
        dge <- calcNormFactors(d, method=norm.method)
    }
    if (verbose) {
        cat("done\n")
    }
    
    # counts per million
    if (verbose) {
        cat("\tGenerating CPM...")
    }
    reads.cpm <- cpm(dge, normalized.lib.sizes=FALSE, log=FALSE, prior.count=0.01)
    # reads.cpm <- cpm(counts.raw, normalized.lib.sizes=FALSE, log=FALSE, prior.count=0.01)
    if (verbose) {
        cat("done\n")
    }
    
    if (export) {
        out.dir <- paste(in.dir,"/bam2cpm-out", sep="")
        if (verbose) {
            cat("\tExporting results to", out.dir, "...")
        }
        dir.create(out.dir, showWarnings = FALSE)
        write.table(cbind(GID=rownames(counts.raw), counts.raw), file=paste(out.dir,"/counts.txt", sep=""), sep="\t", col.names=T, row.names=F)
        write.table(cbind(GID=rownames(reads.cpm), reads.cpm), file=paste(out.dir,"/cpm.txt",sep=""), sep="\t", col.names=T, row.names=F)
        if (verbose) {
            cat("done\n")
        }
    }
    
    return(list(txdb=txdb, ebg=ebg, counts=counts.raw, cpm=reads.cpm))
}


wd <- getwd()


dataset <- "e16.pool1"
in.dir <- paste(wd, "/", dataset, "/", sep="")
print(Sys.time())
ret <- quantifyFromBam(in.dir, feature.type="gene", singleEnd=TRUE, ignore.strand=TRUE, norm.method="none", percentage=1)
print(Sys.time())


dataset <- "e16.pool2"
in.dir <- paste(wd, "/", dataset, "/", sep="")
print(Sys.time())
ret <- quantifyFromBam(in.dir, feature.type="gene", singleEnd=TRUE, ignore.strand=TRUE, norm.method="none", percentage=1)
print(Sys.time())
