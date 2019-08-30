source('~/group/projects/vision/createGenomeTracks.R')

#tss_eRPs_4G = read.table('~/group/projects/vision/rna/tss_eRP_2kb.with0.txt', header=T)
#dist_eRPs_4G = read.table('~/group/projects/vision/rna/dist_eRP_2kb.with0.txt', header=T)
tss_eRPs_4G = read.table('~/group/projects/vision/rna/tss_eRP_2kb.with0.lasso.txt', header=T)
dist_eRPs_4G = read.table('~/group/projects/vision/rna/dist_eRP_2kb.with0.lasso.txt', header=T)


#tss_eRPs_1G = read.table('~/group/projects/vision/rna/tss_eRP_2kb.nosplit4G.with0.txt', header=T)
#dist_eRPs_1G = read.table('~/group/projects/vision/rna/dist_eRP_2kb.nosplit4G.with0.txt', header=T)

tss_eRPs_1G = read.table('~/group/projects/vision/rna/tss_eRP_2kb.nosplit4G.lasso.with0.txt', header=T)
dist_eRPs_1G = read.table('~/group/projects/vision/rna/dist_eRP_2kb.nosplit4G.lasso.with0.txt', header=T)

eRPs0 = read.table('~/group/projects/vision/rna/get_ccRE_gene_pair_nofilter/eRP_all.txt', header=F)


std_range = function(x){
	sx = (x)/(max(x)-min(x))*10
	return(sx)
}

state0_0 = function(x){
	sx = x-x[1]
	return(sx)
}


eRP_all = cbind(eRPs0[,-1], tss_eRPs_1G, dist_eRPs_1G, tss_eRPs_4G, dist_eRPs_4G)

eRP_all = apply(eRP_all, 2, state0_0)

all_rowsum = apply(eRP_all, 1, function(x) sum(x) )
all_colsum = apply(eRP_all, 2, function(x) sum(x) )



dist <- dist(eRP_all, method = "euclidean") # distance matrix
fit <- hclust(dist, method="ward.D2") 
hclust_order = fit$order

library(pheatmap)


colnames(eRP_all) = c('eRP0', 'P0', 'D0', 'P1', 'P2', 'P3', 'P4', 'D1', 'D2', 'D3', 'D4')
rownames(eRP_all) = c(0:26)

corlo_lim_scale = max(max(abs(eRP_all)))
#corlo_lim_scale = 30
breaksList_scale = seq(-corlo_lim_scale, corlo_lim_scale, by = 0.1)
my_colorbar_scale=colorRampPalette(c('blue', 'white', 'red'))(n = length(breaksList_scale))



pdf('eRPs.2kb.all.with0.lasso.pdf', width=4)
pheatmap(eRP_all[order(all_rowsum),], color=my_colorbar_scale, breaks = breaksList_scale, cluster_rows = FALSE, cluster_cols = F,annotation_names_row=FALSE,annotation_names_col=FALSE,show_rownames=T,show_colnames=T)
dev.off()

png('eRPs.2kb.with0.lasso.png', width=200)
pheatmap(eRP_all[order(all_rowsum),], color=my_colorbar_scale, breaks = breaksList_scale, cluster_rows = FALSE, cluster_cols = F,annotation_names_row=FALSE,annotation_names_col=FALSE,show_rownames=T,show_colnames=T)
dev.off()

pdf('eRPs.2kb.hclust.with0.lasso.pdf', width=4)
pheatmap(eRP_all[rev(hclust_order),], color=my_colorbar_scale, breaks = breaksList_scale, cluster_rows = FALSE, cluster_cols = F,annotation_names_row=FALSE,annotation_names_col=FALSE,show_rownames=T,show_colnames=T)
dev.off()




createHeatmap_sort_eRP = function(parafile, statecolor = NULL, markcolor = NULL, cols=c("white","dark blue"), show=TRUE,fout=NULL, sortstate=TRUE,scale=FALSE)
{	x=read.table(parafile,comment="!",header=T);
	k=dim(x)[2];
	l=dim(x)[1];
	p=(sqrt(9+8*(k-1))-3)/2;
	m=as.matrix(x[,1+1:p]/x[,1]);
	colnames(m) = colnames(x)[1+1:p];
	marks=colnames(m);
	rownames(m)=paste(1:l-1," (",round(x[,1]/sum(x[,1])*10000)/100,"%)",sep="");

if(sortstate)
{
o=hclust(dist(m),method="ward.D2")$order;
o=order(-all_rowsum)
m=m[order(-all_rowsum),];
if(length(statecolor) != 0)
{	statecolor=statecolor[o,];
}
}
om=m;
if(scale)
{	m = t((t(m) - apply(m,2,min))/(apply(m,2,max)-apply(m,2,min)+1e-10));
}

	if(length(fout)!=0)
	{	pdf(fout);	}	
	par(mar=c(6,1,1,6));
	rg=range(m);
	colors=0:100/100*(rg[2]-rg[1])+rg[1];
        my_palette=colorRampPalette(cols)(n=100);
	defpalette=palette(my_palette);

if(show)
{
	plot(NA,NA,xlim=c(0,p+0.7),ylim=c(0,l),xaxt="n",yaxt="n",xlab=NA,ylab=NA,frame.plot=F);
	axis(1,at=1:p-0.5,labels=colnames(m),las=2);
	axis(4,at=1:l-0.5,labels=rownames(m),las=2);
	rect(rep(1:p-1,l),rep(1:l-1,each=p),rep(1:p,l),rep(1:l,each=p),col=round((t(m)-rg[1])/(rg[2]-rg[1])*100));#,border=NA);
}
if(scale)
{	m = om;
}

	if(length(statecolor)==0)
	{	if(length(markcolor)==0)
		{	markcolor=t(col2rgb(terrain.colors(ceiling(p))[1:p]));	
			for(i in 1:length(marks))
			{	if(regexpr("h3k4me3",tolower(marks[i]))>0)
				{	markcolor[i,]=c(255,0,0);	}
				if(regexpr("h3k4me2",tolower(marks[i]))>0)
				{	markcolor[i,]=c(250,100,0);	}
				if(regexpr("h3k4me1",tolower(marks[i]))>0)
				{	markcolor[i,]=c(250,250,0);	}
				if(regexpr("h3k36me3",tolower(marks[i]))>0)
				{	markcolor[i,]=c(0,150,0);	}
				if(regexpr("h2a",tolower(marks[i]))>0)
				{	markcolor[i,]=c(0,150,150);	}
				if(regexpr("dnase",tolower(marks[i]))>0)
				{	markcolor[i,]=c(0,200,200);	}
				if(regexpr("atac",tolower(marks[i]))>0)
				{	markcolor[i,]=c(200,50,150);	}
				if(regexpr("h3k9ac",tolower(marks[i]))>0)
				{	markcolor[i,]=c(250,0,200);	}
				if(regexpr("h3k9me3",tolower(marks[i]))>0)
				{	markcolor[i,]=c(100,100,100);	}
				if(regexpr("h3k27ac",tolower(marks[i]))>0)
				{	markcolor[i,]=c(250,150,0);	}
				if(regexpr("h3k27me3",tolower(marks[i]))>0)
				{	markcolor[i,]=c(0,0,225);	}
				if(regexpr("h3k79me2",tolower(marks[i]))>0)
				{	markcolor[i,]=c(200,0,200);	}
				if(regexpr("h4k20me1",tolower(marks[i]))>0)
				{	markcolor[i,]=c(50,200,50);	}
				if(regexpr("ctcf",tolower(marks[i]))>0)
				{	markcolor[i,]=c(200,0,250);	}
			}
		}
		statecolor=array(stateColor(m,markcolor),dim=c(dim(m)[1],2));
	}
	if(show)
	{	rect(rep(p+0.2,l),1:l-0.8,rep(p+0.8,l),1:l-0.2,col=statecolor[,2]);
	}
	if(sortstate)	statecolor[o,]=statecolor;

	palette(defpalette);
	if(length(fout)!=0)
	{	dev.off();	}
	return(statecolor);
}



pdf('pknorm_2_16lim_ref1mo_0424_heatmap.eRP_sort.with0.lasso.pdf', width=4)
createHeatmap_sort_eRP('~/group/projects/vision/pknorm_2_16lim_ref1mo_sample1mo_0424/pknorm_2_16lim_ref1mo_0424_lesshet.para0', sortstate=TRUE)
dev.off()




createHeatmap_sort_eRP = function(parafile, statecolor = NULL, markcolor = NULL, cols=c("white","dark blue"), show=TRUE,fout=NULL, sortstate=TRUE,scale=FALSE)
{	x=read.table(parafile,comment="!",header=T);
	k=dim(x)[2];
	l=dim(x)[1];
	p=(sqrt(9+8*(k-1))-3)/2;
	m=as.matrix(x[,1+1:p]/x[,1]);
	colnames(m) = colnames(x)[1+1:p];
	marks=colnames(m);
	rownames(m)=paste(1:l-1," (",round(x[,1]/sum(x[,1])*10000)/100,"%)",sep="");

if(sortstate)
{
o=hclust(dist(m),method="ward.D2")$order;
o=hclust_order
m=m[hclust_order,];
if(length(statecolor) != 0)
{	statecolor=statecolor[o,];
}
}
om=m;
if(scale)
{	m = t((t(m) - apply(m,2,min))/(apply(m,2,max)-apply(m,2,min)+1e-10));
}

	if(length(fout)!=0)
	{	pdf(fout);	}	
	par(mar=c(6,1,1,6));
	rg=range(m);
	colors=0:100/100*(rg[2]-rg[1])+rg[1];
        my_palette=colorRampPalette(cols)(n=100);
	defpalette=palette(my_palette);

if(show)
{
	plot(NA,NA,xlim=c(0,p+0.7),ylim=c(0,l),xaxt="n",yaxt="n",xlab=NA,ylab=NA,frame.plot=F);
	axis(1,at=1:p-0.5,labels=colnames(m),las=2);
	axis(4,at=1:l-0.5,labels=rownames(m),las=2);
	rect(rep(1:p-1,l),rep(1:l-1,each=p),rep(1:p,l),rep(1:l,each=p),col=round((t(m)-rg[1])/(rg[2]-rg[1])*100));#,border=NA);
}
if(scale)
{	m = om;
}

	if(length(statecolor)==0)
	{	if(length(markcolor)==0)
		{	markcolor=t(col2rgb(terrain.colors(ceiling(p))[1:p]));	
			for(i in 1:length(marks))
			{	if(regexpr("h3k4me3",tolower(marks[i]))>0)
				{	markcolor[i,]=c(255,0,0);	}
				if(regexpr("h3k4me2",tolower(marks[i]))>0)
				{	markcolor[i,]=c(250,100,0);	}
				if(regexpr("h3k4me1",tolower(marks[i]))>0)
				{	markcolor[i,]=c(250,250,0);	}
				if(regexpr("h3k36me3",tolower(marks[i]))>0)
				{	markcolor[i,]=c(0,150,0);	}
				if(regexpr("h2a",tolower(marks[i]))>0)
				{	markcolor[i,]=c(0,150,150);	}
				if(regexpr("dnase",tolower(marks[i]))>0)
				{	markcolor[i,]=c(0,200,200);	}
				if(regexpr("atac",tolower(marks[i]))>0)
				{	markcolor[i,]=c(200,50,150);	}
				if(regexpr("h3k9ac",tolower(marks[i]))>0)
				{	markcolor[i,]=c(250,0,200);	}
				if(regexpr("h3k9me3",tolower(marks[i]))>0)
				{	markcolor[i,]=c(100,100,100);	}
				if(regexpr("h3k27ac",tolower(marks[i]))>0)
				{	markcolor[i,]=c(250,150,0);	}
				if(regexpr("h3k27me3",tolower(marks[i]))>0)
				{	markcolor[i,]=c(0,0,225);	}
				if(regexpr("h3k79me2",tolower(marks[i]))>0)
				{	markcolor[i,]=c(200,0,200);	}
				if(regexpr("h4k20me1",tolower(marks[i]))>0)
				{	markcolor[i,]=c(50,200,50);	}
				if(regexpr("ctcf",tolower(marks[i]))>0)
				{	markcolor[i,]=c(200,0,250);	}
			}
		}
		statecolor=array(stateColor(m,markcolor),dim=c(dim(m)[1],2));
	}
	if(show)
	{	rect(rep(p+0.2,l),1:l-0.8,rep(p+0.8,l),1:l-0.2,col=statecolor[,2]);
	}
	if(sortstate)	statecolor[o,]=statecolor;

	palette(defpalette);
	if(length(fout)!=0)
	{	dev.off();	}
	return(statecolor);
}



pdf('pknorm_2_16lim_ref1mo_0424_heatmap.eRP_sort.hclust.with0.lasso.pdf', width=4)
createHeatmap_sort_eRP('~/group/projects/vision/pknorm_2_16lim_ref1mo_sample1mo_0424/pknorm_2_16lim_ref1mo_0424_lesshet.para0', sortstate=TRUE)
dev.off()










