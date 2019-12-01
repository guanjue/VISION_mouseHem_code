source("degene.R");
#source("genereg.R");
args<-commandArgs(trailingOnly=TRUE);
#method=0;
#if(as.integer(regexpr("chrhmm",args[1]))>0) method=1;
#if(as.integer(regexpr("segway",args[1]))>0) method=2;
#runar2 = TRUE;
#if(runar2) { rungeneReg(args[1], args[2], method=method, mydist=100000); }
#rungeneReg_splinenew(args[1], args[2], method=method);
#runTssEnhEnrich(method,args[1]);

maxdist=1000000;
drange=rbind(c(-1e6,-5e3),c(5e3,1e6));

#direct="dec"
#mypref="step1";
#mymulti=FALSE;
direct="inc"
mypref="step_tt_hp";
mymulti=FALSE;

if(TRUE)
{
DEreg(args[1], as.integer(args[2]), maxdist, mypref,direct,multi = mymulti);
extractDistLoci(paste("/gpfs/scratch/yzz2/rna/",mypref,".",args[2],".",args[1],sep=""));
sampleNull(paste("/gpfs/scratch/yzz2/rna/",mypref,".",args[2],".",args[1],sep=""),args[1],as.integer(args[2]),maxdist=maxdist,fpref=paste(mypref,"null",sep=""));

a=TFenrich(paste("/gpfs/scratch/yzz2/rna/",mypref,".",args[2],".",args[1],sep=""));
write.table(a,paste("/gpfs/scratch/yzz2/rna/",mypref,".",args[2],".",args[1],".tfenrich",sep=""),quote=F,col.names=F);
a=TFenrich(paste("/gpfs/scratch/yzz2/rna/",mypref,".",args[2],".",args[1],sep=""),rmtss=T);
write.table(a,paste("/gpfs/scratch/yzz2/rna/",mypref,".",args[2],".",args[1],".tfenrich_notss",sep=""),quote=F,col.names=F);

a=TFenrich(paste("/gpfs/scratch/yzz2/rna/",mypref,".",args[2],".",args[1],sep=""),factorfile="tsstts.bed");
write.table(a,paste("/gpfs/scratch/yzz2/rna/",mypref,".",args[2],".",args[1],".tssenrich",sep=""),quote=F,col.names=F);
a=TFenrich(paste("/gpfs/scratch/yzz2/rna/",mypref,".",args[2],".",args[1],sep=""),rmtss=T,factorfile="tsstts.bed");
write.table(a,paste("/gpfs/scratch/yzz2/rna/",mypref,".",args[2],".",args[1],".tssenrich_notss",sep=""),quote=F,col.names=F);

a=TFenrich(paste("/gpfs/scratch/yzz2/rna/",mypref,".",args[2],".",args[1],sep=""),distance=250,factorfile="gwas_stopgap_0.05.bed");
write.table(a,paste("/gpfs/scratch/yzz2/rna/",mypref,".",args[2],".",args[1],".gwasenrich",sep=""),quote=F,col.names=F);
a=TFenrich(paste("/gpfs/scratch/yzz2/rna/",mypref,".",args[2],".",args[1],sep=""),distance=250,rmtss=T,factorfile="gwas_stopgap_0.05.bed");
write.table(a,paste("/gpfs/scratch/yzz2/rna/",mypref,".",args[2],".",args[1],".gwasenrich_notss",sep=""),quote=F,col.names=F);


a=eQTLs(paste("/gpfs/scratch/yzz2/rna/",mypref,".",args[2],".",args[1],sep=""),drange=drange);
write.table(a,paste("/gpfs/scratch/yzz2/rna/",mypref,".",args[2],".",args[1],".eqtlenrich",sep=""),quote=F,col.names=F);
a=eQTLs(paste("/gpfs/scratch/yzz2/rna/",mypref,".",args[2],".",args[1],sep=""),drange=drange,rmtss=T);
write.table(a,paste("/gpfs/scratch/yzz2/rna/",mypref,".",args[2],".",args[1],".eqtlenrich_notss",sep=""),quote=F,col.names=F);
}

str0=paste("/gpfs/scratch/yzz2/rna/",mypref,"null.",args[2],".",args[1],".bed",sep="");
if(TRUE)
{
a=TFenrich(paste("/gpfs/scratch/yzz2/rna/",mypref,".",args[2],".",args[1],sep=""),fnull=str0);
write.table(a,paste("/gpfs/scratch/yzz2/rna/",mypref,".",args[2],".",args[1],".tfenrich0",sep=""),quote=F,col.names=F);
a=TFenrich(paste("/gpfs/scratch/yzz2/rna/",mypref,".",args[2],".",args[1],sep=""),rmtss=T,fnull=str0);
write.table(a,paste("/gpfs/scratch/yzz2/rna/",mypref,".",args[2],".",args[1],".tfenrich_notss0",sep=""),quote=F,col.names=F);

a=TFenrich(paste("/gpfs/scratch/yzz2/rna/",mypref,".",args[2],".",args[1],sep=""),fnull=str0,factorfile="tsstts.bed");
write.table(a,paste("/gpfs/scratch/yzz2/rna/",mypref,".",args[2],".",args[1],".tssenrich0",sep=""),quote=F,col.names=F);
a=TFenrich(paste("/gpfs/scratch/yzz2/rna/",mypref,".",args[2],".",args[1],sep=""),rmtss=T,fnull=str0,factorfile="tsstts.bed");
write.table(a,paste("/gpfs/scratch/yzz2/rna/",mypref,".",args[2],".",args[1],".tssenrich_notss0",sep=""),quote=F,col.names=F);

a=TFenrich(paste("/gpfs/scratch/yzz2/rna/",mypref,".",args[2],".",args[1],sep=""),distance=250,fnull=str0,factorfile="gwas_stopgap_0.05.bed");
write.table(a,paste("/gpfs/scratch/yzz2/rna/",mypref,".",args[2],".",args[1],".gwasenrich0",sep=""),quote=F,col.names=F);
a=TFenrich(paste("/gpfs/scratch/yzz2/rna/",mypref,".",args[2],".",args[1],sep=""),rmtss=T,distance=250,fnull=str0,factorfile="gwas_stopgap_0.05.bed");
write.table(a,paste("/gpfs/scratch/yzz2/rna/",mypref,".",args[2],".",args[1],".gwasenrich_notss0",sep=""),quote=F,col.names=F);


a=eQTLs(paste("/gpfs/scratch/yzz2/rna/",mypref,".",args[2],".",args[1],sep=""),drange=drange,fnull=str0);
write.table(a,paste("/gpfs/scratch/yzz2/rna/",mypref,".",args[2],".",args[1],".eqtlenrich0",sep=""),quote=F,col.names=F);
a=eQTLs(paste("/gpfs/scratch/yzz2/rna/",mypref,".",args[2],".",args[1],sep=""),drange=drange,rmtss=T,fnull=str0);
write.table(a,paste("/gpfs/scratch/yzz2/rna/",mypref,".",args[2],".",args[1],".eqtlenrich_notss0",sep=""),quote=F,col.names=F);
}

check(args[1], as.integer(args[2]), fpref=mypref);
checkeqtl(paste("/gpfs/scratch/yzz2/rna/",mypref,".",args[2],".",args[1],sep=""));
check(args[1], as.integer(args[2]), fpref=mypref, pchic = TRUE);
checkeqtl(paste("/gpfs/scratch/yzz2/rna/",mypref,".",args[2],".",args[1],sep=""), suff="pchic");

checkOverlap(paste("/gpfs/scratch/yzz2/rna/",mypref,".",args[2],".",args[1],sep=""), fnull=NULL, eqtl=FALSE, cutoff=5);
checkOverlap(paste("/gpfs/scratch/yzz2/rna/",mypref,".",args[2],".",args[1],sep=""), fnull=paste("/gpfs/scratch/yzz2/rna/",mypref,"null.",args[2],".",args[1],".bed",sep=""), eqtl=FALSE, cutoff=5);
checkOverlap(paste("/gpfs/scratch/yzz2/rna/",mypref,".",args[2],".",args[1],sep=""), fnull=NULL, eqtl=TRUE, cutoff=5);
checkOverlap(paste("/gpfs/scratch/yzz2/rna/",mypref,".",args[2],".",args[1],sep=""), fnull=paste("/gpfs/scratch/yzz2/rna/",mypref,"null.",args[2],".",args[1],".bed",sep=""), eqtl=TRUE, cutoff=5);
#DEpredict(args[1], as.integer(args[2]), maxdist, mypref);
