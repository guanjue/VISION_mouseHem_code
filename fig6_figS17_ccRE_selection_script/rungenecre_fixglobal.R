source("genecre.R");
args<-commandArgs(trailingOnly=TRUE);
rna=read.table("/gpfs/group/yzz2/default/legacy/group/projects/vision/rna/rnaTPM.txt",header=T);
t=which(rna[,1] == args[1]);
x=as.matrix(rna[t,5:16]);
tss=rna[t,2];
m=apply(x,1,mean);
s=apply(x,1,sd);
tt=which(m>-4 & s>2);
e=as.matrix(read.table("vision_genecre.eRP")[,2]);
for(i in 12)
{	rt=run(x[tt,],tss[tt], 12, "../legacy/group/projects/vision/pknorm_2_16lim_ref1mo_sample1mo_0424/pknorm_2_16lim_ref1mo_0424_lesshet",args[1],e,lessone=i,B=100,fixEffect=TRUE);
	save(rt,file=paste("vision_result/vision_rna_tss5k_fix_global.",args[1],".3.",i,".Rdat",sep=""));
}
tt=which(m>-4 & s<=2);
for(i in 1:12)
{	rt=run(x[tt,],tss[tt], 12, "../legacy/group/projects/vision/pknorm_2_16lim_ref1mo_sample1mo_0424/pknorm_2_16lim_ref1mo_0424_lesshet",args[1],e,lessone=i,B=100,fixEffect=TRUE);
	save(rt,file=paste("vision_result/vision_rna_tss5k_fix_global.",args[1],".4.",i,".Rdat",sep=""));
}
tt=which(m<=-4 & s>2);
for(i in 1:12)
{	rt=run(x[tt,],tss[tt], 12, "../legacy/group/projects/vision/pknorm_2_16lim_ref1mo_sample1mo_0424/pknorm_2_16lim_ref1mo_0424_lesshet",args[1],e,lessone=i,B=100,fixEffect=TRUE);
	save(rt,file=paste("vision_result/vision_rna_tss5k_fix_global.",args[1],".2.",i,".Rdat",sep=""));
}
tt=which(m<=-4 & s<=2);
for(i in 1:12)
{	rt=run(x[tt,],tss[tt], 12, "../legacy/group/projects/vision/pknorm_2_16lim_ref1mo_sample1mo_0424/pknorm_2_16lim_ref1mo_0424_lesshet",args[1],e,lessone=i,B=100,fixEffect=TRUE);
	save(rt,file=paste("vision_result/vision_rna_tss5k_fix_global.",args[1],".1.",i,".Rdat",sep=""));
}
