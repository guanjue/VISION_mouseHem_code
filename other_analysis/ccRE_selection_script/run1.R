source("degene_try.R");
args<-commandArgs(trailingOnly=TRUE);
#DEpredict(args[1], "/gpfs/group/yzz2/default/scratch/roadmap_analysis/impute/","bin98.",fhic="../../legacy/group/HiC/Gm12878/10kb/chr16.matrix",hicresolution=1e4,fpref="tt",maxdist=100000,gidlist=as.integer(args[2]):as.integer(args[3]));
#DEpredict(args[1], "/gpfs/group/yzz2/default/scratch/roadmap_analysis/impute/","bin98.",fhic="../../legacy/group/HiC/Gm12878/10kb/chr16.matrix",hicresolution=1e4,fpref="tt1m",maxdist=1000000,gidlist=as.integer(args[2]):as.integer(args[3]));
DEpredict(args[1], "/gpfs/group/yzz2/default/scratch/roadmap_analysis/impute/","bin98.",fpref="n1m",maxdist=1000000,gidlist=as.integer(args[2]):as.integer(args[3]));
