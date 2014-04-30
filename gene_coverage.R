#!/usr/bin/Rscript
# Rscript for coverage statistics for a list of genes in a sample
# Rscript coverage.R '<full path to bam file>' geneA geneB geneX minCoverage


VERSION <- '1.3-Mar2014'

#suppressMessages(require(multicore,quiet=TRUE))
#suppressMessages(require(GenomicFeatures,quiet=TRUE))

# get the curr dir and source shared functions
args <- commandArgs(trailingOnly = FALSE)
script.basename <- dirname(sub('--file=', '', args[grep('--file=', args)]))
source(paste(script.basename, 'shared_functions.R',sep='/'))

args <- commandArgs(trailingOnly = TRUE)
bam <- args[1]
genes <- args[2:length(args)]

cat('Script version', VERSION,'\n')
cat('Sample', bam,'\n')
cat('Candidate Genes: ')
cat(genes)
cat('\n\n')
cat("Individual gene coverage analysis")
cat('\n')

#iterate over each gene
for (gene in genes) {
 
  stats <- get.gene.coverage.stats(gene,bam)
    
  if (dim(stats[['ensembl_gene_record']])[1] > 0) {

    cvrg.df <- stats[['coverage.stats']]
     
    cat('Gene ', gene," - ",stats[['number.transcripts']]," transcripts, ",stats[['number.exons']]," exons, ", 
        formatC(cvrg.df[['entire gene']][1],digits=2,format='f'),'X mean exonic coverage, ', 
        formatC(cvrg.df[['entire gene']][2],digits=2,format='f'),' above 10X, ', 
        formatC(cvrg.df[['entire gene']][3],digits=2,format='f'),' above 20X, ', 
        formatC(cvrg.df[['entire gene']][4],digits=2,format='f'),' above 30X',
        sep='')
    cat("\n")
    cat("List of exons:")
    cat("\n")
    cat(paste('chromosome','start','end', 'mean_coverage','%bp>10X','%bp>20X','%bp>30X','\n',sep='\t'))
    
    # coverage is calculated for unique exon cooridnates
    exon.coords<-unique(stats[['ensembl_gene_record']][,c('chromosome_name','genomic_coding_start','genomic_coding_end')])
    
    for (j in 2:ncol(cvrg.df)) {
      cat(paste(exon.coords$chromosome_name[1], 
                exon.coords$genomic_coding_start[j-1], 
                exon.coords$genomic_coding_end[j-1], 
                round(cvrg.df[1,j],2),
                round(cvrg.df[2,j],2)*100, round(cvrg.df[3,j],2)*100, round(cvrg.df[4,j],2)*100, '\n', sep="\t"))
    }
    
    cat('\n')
  } else {
    cat(paste("Gene symbol",gene,"not found."))
    cat('\n')
  }
}

