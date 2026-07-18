#!/usr/bin/env Rscript

## run example:
##  Rscript Downsample_fragmentfile.R /proj/yunligrp/users/lagler/HiC/fithic/GM12878_10Kb/downsampled/10percent/ MAPQGE30_chr _fhic2.downsampled_0.1.txt.gz _fhic1.downsampled_0.1
##
## arguments:
## INFDIR - dir with reg files
## PREFIX - dataset prefix (reg_raw for example)
## SUFFIX - dataset suffix
## chroms - number of shromosomes (19 for mouse, 22 for human)
## OUTSUFIX - suffix for output file name

library(data.table)
library(R.utils)
args = commandArgs(trailingOnly = T)

if (length(args) != 5) {
  print('Wrong number of arguments. Stopping.')
  print('Arguments needed (in this order): INFDIR, PREFIX, SUFFIX, chroms, OUTSUFFIX')
  print(paste('Number of arguments entered:',length(args)))
  print('Arguments entered:')
  print(args)
  quit()
} else {
  print(args)
  INFDIR = args[1]
  PREFIX = args[2]
  SUFFIX = args[3]
  chroms = args[4]
  OUTSUFFIX = args[5]
}

options(scipen=999)

# Load data, sample counts, split data, write downsampled datasets
for (i in 1:as.numeric(chroms)) {
  print(paste0('loading chromosome ',i))
  mm <- fread(paste0(INFDIR,PREFIX,i,SUFFIX), header = F, stringsAsFactors = F, data.table = T)
  colnames(mm) <- c("chr1", "fragment", "ch2", "frag2", "count")
  
  mm2 <- mm[,list(mcc = sum(count)), by = c("fragment")] #marginalized contact count by fragment mid

  n.frag <- length(mm2$fragment)
  chr <- as.vector(rep(i, n.frag)) #chromosome
  
  map <- NULL; ef <- NULL; 
  ef <- as.vector(rep(0, n.frag)) #extraField
  map = as.integer(mm2$mcc >= 1) #mappable?
  
  fhic1 <- cbind(chr,ef, mm2$fragment, mm2$mcc, map)
  fhic1_df <- data.frame(fhic1)

  print(paste0('Writing downsampled dataset ',i))
  fwrite(fhic1, paste0(INFDIR,PREFIX,i,OUTSUFFIX, ".txt"), sep="\t", row.names = F, col.names = F, quote=F)
  gzip(paste0(INFDIR,PREFIX,i,OUTSUFFIX,".txt"),
       destname=paste0(INFDIR,PREFIX,i,OUTSUFFIX,".txt.gz"),
       overwrite=TRUE)
}



