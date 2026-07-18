#!/usr/bin/env Rscript

## run example:
##  Rscript Downsample.R /proj/yunligrp/users/lagler/HiC/fithic/GM12878_10Kb/FitHiC_inputs/ MAPQGE30_ _raw_count_matrix.txt.fhic2.gz 22 0.1 /proj/yunligrp/users/lagler/HiC/fithic/GM12878_10Kb/downsampled/10percent/
##
## arguments:
## INFDIR - dir with reg files
## PREFIX - dataset prefix (reg_raw for example)
## SUFFIX - dataset suffix
## chroms - number of shromosomes (19 for mouse, 22 for human)
## DRatio - Downsampling depth ratio relative to original depth
## OUTDIR - output directory

library(data.table)
library(R.utils)
args = commandArgs(trailingOnly = T)

if (length(args) != 7) {
  print('Wrong number of arguments. Stopping.')
  print('Arguments needed (in this order): INFDIR, PREFIX, SUFFIX, chroms, DRatio, DRatio_filename, OUTDIR')
  print(paste('Number of arguments entered:',length(args)))
  print('Arguments entered:')
  print(args)
  quit()
} else {
  print(args)
  INFDIR = args[1]
  PREFIX = args[2]
  SUFFIX = args[3]
  chroms = paste0('chr',seq(1,as.numeric(args[4]),1))
  DRatio = args[5]
  DRatio_filename = args[6]
  OUTDIR = args[7]
}


# Load data, sample counts, split data, write downsampled datasets
for (i in chroms) {
  print(paste0('loading chromosome ',i))
  mm = fread(paste0(INFDIR,PREFIX,i,SUFFIX), header = F, stringsAsFactors = F, data.table = F)
  n.count = length(mm[,5])
  n.reads = sum(mm[,5])
  mm[,5] = rmultinom(1, floor(as.numeric(DRatio) * n.reads), mm[,5] / n.reads)	
  print(paste0('Writing downsampled dataset ',i))
  fwrite(mm, paste0(OUTDIR,PREFIX,i,'_fhic2.downsampled_', DRatio_filename, '.txt'), sep="\t", row.names = F, col.names = F, quote=F)
  gzip(paste0(OUTDIR, PREFIX,i,'_fhic2.downsampled_', DRatio_filename, '.txt'),
       destname=paste0(OUTDIR,PREFIX,i,'_fhic2.downsampled_', DRatio_filename, '.txt.gz'),
       overwrite=TRUE)
}



