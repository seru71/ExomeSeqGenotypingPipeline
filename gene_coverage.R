#!/usr/bin/Rscript
# Rscript for coverage statistics for a list of genes in a sample
# Rscript coverage.R '<full path to bam file>' geneA geneB geneX minCoverage

VERSION <- '1.3-Mar2014'

.libPaths('/export/astrakanfs/stefanj/R/library')
suppressMessages(require(Rsamtools,quiet=TRUE))
suppressMessages(require(multicore,quiet=TRUE))
suppressMessages(require(GenomicFeatures,quiet=TRUE))
suppressMessages(require(GenomicRanges,quiet=TRUE))
suppressMessages(require(biomaRt,quiet=TRUE))
options(error=traceback)

# get the curr dir and source shared functions
args <- commandArgs(trailingOnly = FALSE)
script.basename <- dirname(sub('--file=', '', args[grep('--file=', args)]))
source(paste(script.basename, 'shared_functions.R',sep='/'))

args <- commandArgs(trailingOnly = TRUE)
bam <- args[1]
genes <- args[2:length(args)]

# create index file if it does not exist
createBamIndex(bam)

#load ensembl mart
human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
genes_info <- getBM(
      attributes = c("external_gene_id", "ensembl_gene_id", "ensembl_transcript_id", 
                     "chromosome_name", "genomic_coding_start","genomic_coding_end", "strand"), 
      filters = c("hgnc_symbol","with_ccds"), 
      values = list(hgnc_symbol = genes, with_ccds=TRUE), 
      mart=human)

cat('Script version', VERSION,'\n')
cat('Sample', bam,'\n')
cat('Candidate Genes: ')
cat(genes)
cat('\n\n')
cat("Individual gene coverage analysis")
cat('\n')

#iterate over each gene
for (i in 1:length(genes)) {
  info <- subset(genes_info, external_gene_id == genes[i])
  # Remove any references to Locus Reference Genomic locations
  info <- subset(info, !grepl("LRG", ensembl_gene_id))
  
  if (dim(info)[1] > 0) {
    #focus on the gene region
    chr <- as.character(info$chromosome_name[1]) # assuming that all gene exons are on one chromosome
    bamRegion <- getSpecificRegion(chr, min(info$genomic_coding_start), max(info$genomic_coding_end), bam)
    
    exons <- GRanges(chr,IRanges(info$genomic_coding_start, info$genomic_coding_end), strand=info$strand)
    # remove identical exons
    exons <- unique(exons)
    #seqlevels(exons) <- sub("^(\\d+)","chr\\1",seqlevels(exons))
    
    number.transcripts <- length(unique(info$ensembl_transcript_id))   
    # incorrect, should calc unique ranges in exons object
    number.exons <- nrow(as.data.frame(exons))
    
    #intersect the transcript range with the actual reads reported, to calculate coverage
    coverage.exons <- Views(coverage(bamRegion)[chr], as(exons,"RangesList")[chr])[[1]] #, width=max(info$genomic_coding_end)
    
    # calc coverage on merged overlapping exons for the purpose of correct gene-wide cvrg calculation
    coverage.merged.exons <- Views(coverage(bamRegion)[chr], as(reduce(exons),"RangesList")[chr])[[1]]
    
    #The mean is the weighted mean
    gene.mean.coverage <- viewMeans(coverage.merged.exons)%*%width(coverage.merged.exons)/sum(width(coverage.merged.exons))

    gene.length = sum(width(coverage.merged.exons))
    X10 <- sum(viewApply(coverage.merged.exons, function(x) sum(width(slice(x,lower=10))))) / gene.length
    X20 <- sum(viewApply(coverage.merged.exons, function(x) sum(width(slice(x,lower=20))))) / gene.length
    X30 <- sum(viewApply(coverage.merged.exons, function(x) sum(width(slice(x,lower=30))))) / gene.length
    cat('Gene ', genes[i]," - ",number.transcripts," transcripts, ",number.exons," exons, ", 
        formatC(gene.mean.coverage,digits=2,format='f'),'X mean exonic coverage, ', 
        formatC(X10,digits=2,format='f'),' above 10X, ', 
        formatC(X20,digits=2,format='f'),' above 20X, ', 
        formatC(X30,digits=2,format='f'),' above 30X',
        sep='')
    cat("\n")
    cat("List of exons:")
    cat("\n")
    cat(paste('chromosome','start','end', 'mean_coverage','>10X','>20X','>30X','\n',sep='\t'))
    
    for (j in 1:number.exons) {
      meanCoverage <- formatC(as.vector(viewMeans(coverage.exons))[j],digits=2,format="f")
      exon.start <- as.data.frame(exons)[j,][2]
      exon.end <- as.data.frame(exons)[j,][3]
      exon.size <- sum(width(coverage.exons[j]))
      exon.10X <- sum(width(slice(coverage.exons[[j]], lower=10))) / exon.size
      exon.20X <- sum(width(slice(coverage.exons[[j]], lower=20))) / exon.size
      exon.30X <- sum(width(slice(coverage.exons[[j]], lower=30))) / exon.size
      cat(paste(chr, exon.start,exon.end,meanCoverage,round(exon.10X, 2),round(exon.20X, 2) ,round(exon.30X, 2), '\n', sep="\t"))
    }
    
    cat('\n')
  } else {
    cat(paste("Gene symbol",genes[i],"not found."))
    cat('\n')
  }
}
