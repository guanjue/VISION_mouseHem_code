
source('../../createGenomeTracks.R')



createHeatmap<-function(parafile, statecolor = NULL, markcolor = NULL, cols=c("white","dark blue"), show=TRUE,fout=NULL, sortstate=TRUE,scale=FALSE)
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
m=m[o,];
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
#if(scale)
#{	m = om;
#}
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
				if(regexpr("dnase",tolower(marks[i]))>0)
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
				if(regexpr("wgbs",tolower(marks[i]))>0)
				{	markcolor[i,]=c(30,144,255);	}
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


pdf('run_IDEAS.para.check.pdf')
sc=createHeatmap('run_IDEAS.para.check',scale=F)
dev.off()

pdf('vision_withMeth.pdf')
sc=createHeatmap('run_IDEAS.para',scale=F)
dev.off()


createTrack('run_IDEAS.reordered.state', 'mm10.genome', 'run_IDEAS', sc, NULL, NULL)

