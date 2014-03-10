

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


  # if bam has chr coordinates with "chr" prefix
  if (length(grep("chr",seqnames(seqinfo(BamFile(bamFile))),ignore.case=T)) > 0) {    
    chr <- paste("chr",chr,sep="") 
  }

  param <- ScanBamParam(what = c("rname", "strand","pos", "qwidth"),
                        which = GRanges(chr,IRanges(chrStart, chrEnd)),
                        flag = scanBamFlag(isUnmappedQuery = FALSE)
  )                    
  
  x <- scanBam(bamFile, param = param)[[1]]
  ranges = GRanges(seqnames=Rle(x$rname), ranges=IRanges(x$pos, width=x$qwidth))
  seqlevels(ranges) <- sub("^(\\d+)","chr\\1",seqlevels(ranges))
  ranges
  # coverage(ranges)
  #coverage(IRanges(x[["pos"]], width = x[["qwidth"]]))
  ### coverage IRLe object
}

