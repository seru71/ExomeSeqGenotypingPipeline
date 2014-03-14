
if (file.exists('/export/astrakanfs')) {
  .libPaths('/export/astrakanfs/stefanj/R/library')
}
suppressMessages(require(Rsamtools,quiet=TRUE))
suppressMessages(require(GenomicRanges,quiet=TRUE))
suppressMessages(require(biomaRt,quiet=TRUE))
options(error=traceback)


#
# creates index for bam file if it does not exist
createBamIndex <- function(bamFile  ### path to the bam file
                          ) {
  baiFile = sub(".bam$",".bai",bamFile)
  if (!file.exists(baiFile)) {
    indexBam(bamFile)
    file.rename(paste(bamFile, '.bai', sep=''), baiFile)
  } 
}


#
# function to return a GRanges object for specific region of a bam file
# should reduce memory consumption because the bam file is not read entirely
getSpecificRegion <- function(chr, ### chromosome
                              chrStart, ### start position (bp)
                              chrEnd,  ### end position (bp)
                              bamFile ### bam file
                              ) {
  
  createBamIndex(bamFile)

  # if bam has chromosome names with "chr" prefix
  if (length(grep("chr",seqnames(seqinfo(BamFile(bamFile))),ignore.case=T)) > 0) {    
    chr <- paste("chr",chr,sep="") 
  }
  
  param <- ScanBamParam(what = c("rname", "strand","pos", "qwidth"),
                        which = GRanges(chr,IRanges(chrStart, chrEnd)),
                        flag = scanBamFlag(isUnmappedQuery = FALSE)
  )                    
  
  x <- scanBam(bamFile, param = param)[[1]]
  ranges = GRanges(seqnames=Rle(x$rname), ranges=IRanges(x$pos, width=x$qwidth))
  #seqlevels(ranges) <- sub("^(\\d+)","chr\\1",seqlevels(ranges))
  seqlevels(ranges) <- sub("^chr","",seqlevels(ranges))
  return(ranges)
}

# global vars, they should be part of an object!!
human=NA

get.gene.info <- function(gene) {
  
  #load ensembl mart
  if (!isS4(human)) {
    cat('Getting mart...')
    human <<- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
    cat('done.\n')
  }
  gene_info <- getBM(
      attributes = c("external_gene_id", "ensembl_gene_id", "ensembl_transcript_id", 
                     "chromosome_name", "genomic_coding_start","genomic_coding_end", "strand"), 
      filters = c("hgnc_symbol","with_ccds"), 
      values = list(hgnc_symbol = gene, with_ccds=TRUE), 
      mart=human)
  
  info <- subset(gene_info, external_gene_id == gene)
  # Remove any references to Locus Reference Genomic locations
  info <- subset(info, !grepl("LRG", ensembl_gene_id))
  return(info)
}


get.gene.coverage.stats <- function(gene, bam) {
  
  info <- get.gene.info(gene)
  coverage.stats <- list(ensembl_gene_record=info)
  
  if (dim(info)[1] > 0) {
    #focus on the gene region
    chr <- as.character(info$chromosome_name[1]) # assuming that all gene exons are on one chromosome
    bamRegion <- getSpecificRegion(chr, min(info$genomic_coding_start), max(info$genomic_coding_end), bam)
    
    exons <- GRanges(chr,IRanges(info$genomic_coding_start, info$genomic_coding_end), strand=info$strand)
    # remove identical exons
    exons <- unique(exons)
    #seqlevels(exons) <- sub("^(\\d+)","chr\\1",seqlevels(exons))
    
    coverage.stats[['number.transcripts']] <- length(unique(info$ensembl_transcript_id))   
    # incorrect, should calc unique ranges in exons object
    coverage.stats[['number.exons']] <- nrow(as.data.frame(exons))
            

    cvrg.df <- data.frame(row.names=c("mean coverage","% bp above 10X","% bp above 20X","% bp above 30X"))    
    
    #
    # whole gene stats
    #
    # calc coverage on merged overlapping exons for the purpose of correct gene-wide cvrg calculation
    coverage.merged.exons <- Views(coverage(bamRegion)[chr], as(reduce(exons),"RangesList")[chr])[[1]]
    gene.length = sum(width(coverage.merged.exons))
    
    #The gene coverage mean is a weighted mean
    cvrg.df['entire gene'] <- c(viewMeans(coverage.merged.exons)%*%width(coverage.merged.exons)/sum(width(coverage.merged.exons)),
                                sum(viewApply(coverage.merged.exons, function(x) sum(width(slice(x,lower=10))))) / gene.length,
                                sum(viewApply(coverage.merged.exons, function(x) sum(width(slice(x,lower=20))))) / gene.length,
                                sum(viewApply(coverage.merged.exons, function(x) sum(width(slice(x,lower=30))))) / gene.length)
     
    #
    # per exome stats
    #
    #intersect the transcript range with the actual reads reported, to calculate coverage
    coverage.exons <- Views(coverage(bamRegion)[chr], as(exons,"RangesList")[chr])[[1]] #, width=max(info$genomic_coding_end)
    
    for (j in 1:coverage.stats[['number.exons']]) {
      exon.start <- as.data.frame(exons)[j,][2]
      exon.end <- as.data.frame(exons)[j,][3]
      exon.size <- sum(width(coverage.exons[j]))
      
      cvrg.df[paste("exon",j)] <- c(as.vector(viewMeans(coverage.exons))[j],
                                    sum(width(slice(coverage.exons[[j]], lower=10))) / exon.size,
                                    sum(width(slice(coverage.exons[[j]], lower=20))) / exon.size,
                                    sum(width(slice(coverage.exons[[j]], lower=30))) / exon.size)
    }
    
    coverage.stats[['coverage.stats']] <- cvrg.df
#    coverage.stats[['exon.starts']] <- as.data.frame(exons)['start']
#    coverage.stats[['exon.ends']] <- as.data.frame(exons)['end']
    
  }
  
  return(coverage.stats)
}

