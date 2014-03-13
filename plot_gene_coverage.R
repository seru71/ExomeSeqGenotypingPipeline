#!/usr/bin/Rscript
args <- commandArgs(trailingOnly = FALSE)
script.basename <- dirname(sub('--file=', '', args[grep('--file=', args)]))
source(paste(script.basename, 'shared_functions.R',sep='/'))

genes=c("HNF1A","HNF4A")
#bams=c("2351-SJ-0001_S1.bam","2351-SJ-0002_S2.bam","2351-SJ-0003_S3.bam")
bams=paste("/export/astrakanfs/stefanj/data/amplicon_based_Jan2014/",
           c(paste("2351-SJ-000",1:9,"_S",1:9,sep=""), paste("2351-SJ-00",10:96,"_S",10:96,sep="")),
           ".bam",sep="")

mean.gene.cvrg = list()  

for (gene in genes) {
  
  cat(paste("Analyzing ", gene, "...\n", sep=""))
  
  figure.name=paste(gene,'-exon-coverage-stats.png',sep='')
  png(figure.name)
  
  for (bam in bams) {
    
    cat(paste("in", bam,"\n")) 
        
    stats <- get.gene.coverage.stats(gene,bam)[['coverage.stats']]
    
    number.exons = ncol(stats)-1
    boxplot.data = c()
    line.stat = "% bp above 20X"
    
    Y_MAX=1500
    if (bam == bams[1]) {  # the first bam
      boxplot(as.list(rep(NA,number.exons)), mar=c(5,4,4,2), xlab="Exons", ylim=c(1,Y_MAX), ylab="Coverage")    
      axis(4,at=c(0:10)*(Y_MAX/10), labels=c(0:10)*0.1, col='red')
      mtext(line.stat, 4)
    } 
    boxplot.data = c(boxplot.data, unlist(stats['mean coverage',(1:number.exons)+1]))   # mean cvrg for all exons
    points(unlist(stats[line.stat,(1:number.exons)+1]) * Y_MAX, type='l', col='red')  # line.stat statistic for all exons
    
    # for all-genes-boxplot
    mean.gene.cvrg[[gene]] <- c(mean.gene.cvrg[[gene]], stats['mean coverage','entire gene'])
  }
  
  boxplot(matrix(boxplot.data, byrow=T, nc=number.exons), add=T)
  title(gene)
  dev.off()    
  
  cat(paste("Saved ", figure.name))
}

png('all-genes-mean-coverage.png')
boxplot(mean.gene.cvrg)
dev.off()
cat("Saved all-genes-mean-coverage.png")

