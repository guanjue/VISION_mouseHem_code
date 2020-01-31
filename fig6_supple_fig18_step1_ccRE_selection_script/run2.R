source("degene_try.R");
args<-commandArgs(trailingOnly=TRUE);
DEpredict(args[1], "/gpfs/scratch/yzz2/encode_mouse/","me66n.",frna="me66_gene_tpm_full.txt",fpref="me66n1m",maxdist=1000000,gidlist=as.integer(args[2]):as.integer(args[3]));
