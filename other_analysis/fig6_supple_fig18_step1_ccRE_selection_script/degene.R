#DEreg: use stepwise regression to find interacting loci
#extractDistLoci: generate position index of loci and save in a file
#TFenrich: use encode TF binding cluster to perform enrichment analysis
#stepMerge: mcmc version of finding sum regression variables
#stepMerge_inc, stepMerge_dec: greedy stepwise selection for sum regression
#sampleNULL: sample variables randomly with matched counts per gene corresponding to cases

#run sequence: DEreg -> extractDistLoci (-> sampleNull) -> TFenrich

minr = 0.2;
mpenalty=0.0;#0.5;1;

DEreg<-function(chr, path, statepref, usemean = F, maxdist=200000, fpref="step", direction="inc", resume=NULL, multi=FALSE)
{	gene=read.table("Ensembl_v65.Gencode_v10.ENSG.gene_info");
	Y=array(0,dim=c(dim(gene)[1],56));
	x=read.table("57epigenomes.RPKM.pc",header=T);
	cells=colnames(x)[-(1:2)];
	m=match(gene[,1],x[,1]);
	t=which(is.na(m)==F);
	Y[t,]=as.matrix(x[m[t],-(1:2)]);
	x=read.table("57epigenomes.RPKM.nc",header=T);
	m=match(gene[,1],x[,1]);
	t=which(is.na(m)==F);
	Y[t,]=as.matrix(x[m[t],-(1:2)]);
	x=read.table("57epigenomes.RPKM.rb",header=T);
	m=match(gene[,1],x[,1]);
	t=which(is.na(m)==F);
	Y[t,]=as.matrix(x[m[t],-(1:2)]);

	Y=log2(Y+0.01);

	t=which(apply(Y,1,sd) > 0);
	Y=Y[t,];
	gene=gene[t,];	
	genechr = paste("chr",gene[,2],sep="");
	tss = gene[,3];
	tss[gene[,5]<0] = gene[gene[,5]<0,4];

	library("data.table");
	library("MASS");
	cellinfo = read.table("/gpfs/group/yzz2/default/RoadmapEpigenomics/extra/epigenomeID.txt",header=T,comment="!",sep="\t");
	cellinfo = cellinfo[order(as.integer(substr(cellinfo[,1],2,1000))),];
	cellmap = match(cells, cellinfo[,1]);
	nn = dim(Y)[2];
	for(i in chr)
	{	print(i);
		Z = NULL;
		g = fread(paste(path,statepref,i,".state",sep=""));
		a=system(paste("wc -l ",path, statepref,i,".para",sep=""),intern=T);
		gn=as.integer(substr(a,1,regexpr(" ",a)[1]))-1;
		pos = as.matrix(g[,3]);
		state = as.matrix(g[,5:131])[,cellmap];
		t = which(gene[,2]==substr(i,4,1000));

		if(usemean)
		{	p0 = read.table(paste(path,statepref,i,".para",sep=""));
			p0 = as.matrix(p0[,2:13]/p0[,1]);
			Z = t(array(c(t(p0[c(t(state))+1,])),dim=c(length(cellmap)*12,dim(state)[1])));
		}

		print("state file loaded");
	
		gid = tinv = rinv = NULL;	
		for(j in t)
		{	tt = which(pos >= tss[j] - 1000 & pos < tss[j] + 1000);
			if(length(tt) == 0) next;
			rr = which(pos >= tss[j] - maxdist & pos < tss[j] + maxdist);
			if(length(rr) < 2) next;
			gid = c(gid, j);
			tinv = rbind(tinv, range(tt));
			rinv = rbind(rinv, range(rr));	
		}

		if(length(gid) < 10) next;
		ee=NULL;
		str = paste(path,statepref,"coef",sep="");
		if(file.exists(str) == F)
		{	ee = getStateCoef(paste(path,statepref,sep=""));
		} else
		{	ee = read.table(str);
		}
		
		lab = rownames(ee);
		l = apply(t(t(lab)),1,function(x){tail(unlist(gregexpr("_",x)),n=1)});
		lab = substr(lab,1,l-1);
		ee = t(array(ee[,(dim(ee)[2]+1)/2],dim=c(gn,dim(ee)[1]/gn)));
		rownames(ee) = unique(lab);
		colnames(ee) = 1:gn-1;
		t=which(apply(ee,1,sd) == 0);
		if(length(t)>0) ee[t,] = ee[dim(ee)[1],];

		fout = paste("/gpfs/scratch/yzz2/rna/", fpref, ".", statepref,i,sep="");
		gst = 0;
		if(length(resume)>0) gst=which(gene[gid,7]==resume)[1];
		for(j in (gst+1):length(gid))
		{	
			if(length(ee)>0) e = ee[which(gene[gid[j],6]==rownames(ee)),];	
			X = t(array(e[c(state[rinv[j,1]:rinv[j,2],])+1], dim=c(rinv[j,2]-rinv[j,1]+1,nn)));
			#X = smoothMatrix(X);
			yy = Y[gid[j],];
			r = cor(yy, X+rnorm(length(X),0,0.0001));
			r[is.na(r)==T] = 0;
			a = which(r > minr);
#a=which((r-mean(r))/sd(r) > qnorm(1-0.2));
			if(usemean)
			{	bb = NULL;
				for(jj in tinv[j,1]:tinv[j,2])
				{	trr = cor(Z[jj,],t(Z[rinv[j,1]:rinv[j,2],]))[1,];
					trr[is.na(trr)==T]=0;
					trr=(trr-mean(trr))/sd(trr);
					bb=c(bb,which(trr>qnorm(1-0.05/(tinv[j,2]-tinv[j,1]+1))));
				}
				bb = sort(unique(bb));
				bbb=NULL;
				jj = 1; 
				for(kk in 2:length(bb))
				{	if(bb[kk]-bb[kk-1]>10 | kk==length(bb))
					{	bbb=rbind(bbb,c(bb[jj],bb[kk-as.integer(kk<length(bb))]));
						jj = kk;
					}
				}
				bbb=array(bbb,dim=c(length(bbb)/2,2));
				bb=NULL;for(kk in 1:dim(bbb)[1]) bb=sort(unique(c(bb,bbb[kk,1]:bbb[kk,2])));
				#print(c(length(a),length(bb),length(unique(c(a,bb)))));
				#a = sort(unique(c(a, bb)));
				print(c(length(a),length(bb),length(intersect(a,bb))));
				a = sort(intersect(a, bb));
			}
			#r = r^2;

#			a = order(r)[-(1:(length(r)-min(40,length(r)-1)))];
#	a = a[which(r[a]>0.1)];
#			b = tinv[j,1]:tinv[j,2]-rinv[j,1]+1;
#	b=b[which(r[b]>0.1)];
#			a = sort(unique(c(a, b)));
#	if(length(a)==0 | length(b)==0) next;
#			str=paste("yy~",paste(paste("X[,",a,"]",sep=""),collapse="+"),sep="");
#			btr=paste("~",paste(paste("X[,",b,"]",sep=""),collapse="+"),sep="");

#			bmax = b[order(r[b])[length(b)]];
#			rr = stepAIC(lm(paste("yy",btr,sep="")),direction="both",k=log(nn),scope=list(lower= paste("~X[,",bmax,"]",sep="")), trace=F);	
#			if(summary(rr)$adj.r.squared < 0.01) 
#			{	b = bmax; 
#			} else
#			{	b = as.integer(unlist(strsplit(paste(as.character(rr$terms[[3]]),collapse=" "),"[^0-9]+"))); 
#				b = b[is.na(b)==F];
#			}
#			btr=paste("~",paste(paste("X[,",b,"]",sep=""),collapse="+"),sep="");
#
#			rr = stepAIC(lm(str),direction="both",k=log(nn),scope=list(lower= btr),trace=F);
#			r2 = summary(rr)$adj.r.squared;
#			if(r2 < 0.01) next;
#			ii = as.integer(unlist(strsplit(paste(as.character(rr$terms[[3]]),collapse=" "),"[^0-9]+")));
#			ii = ii[is.na(ii)==F];

			#b = tinv[j,1]:tinv[j,2] - rinv[j,1] + 1;
			#ii = stepMerge_inc(array(X[,b],dim=c(nn, length(b))), yy);
			#ii = b[ii];
			#a = sort(unique(c(which(r > 0.1), b)));;
			#ii = stepMerge_inc(array(X[,a],dim=c(nn, length(a))), yy, sel = match(ii, a));
			#ii = a[ii];
			
			#a = 1:dim(X)[2];
			if(length(a) < 1) next;
			#rti = runStepMerge(X[,a], yy, direction = direction);
			if(multi)
			{	if(direction == "inc") {	rti = stepMulti_inc(X[,a], yy, mpenalty);
				} else
				{	rti = stepMulti_dec(X[,a],yy,mpenalty);
				}
				mcmc=cbind(rti,5);
				ss=rti;
				rti=NULL;rti$mcmc=mcmc;rti$ss=ss;	
			} else
			{	rti = runStepMerge1(X[,a], yy, penalty=mpenalty, direction = direction);
			}
			if(length(rti$ss) < 1) next;
			ii = a[rti$ss];
			rti$mcmc[,1] = a[rti$mcmc[,1]];
#			b = tinv[j,1]:tinv[j,2] - rinv[j,1] + 1;
#			if(length(b) < 1) next;
#			if(length(intersect(ii, b)) == 0)
#			{	xsum = apply(array(X[,ii],dim=c(nn,length(ii))),1,sum);
#				olp = -1000000000;
#				iii = b[1];
#				for(jj in b)
#				{	tlp = logLik(lm(yy~I(xsum+X[,jj])));
#					if(jj == b[1] | tlp > olp) 
#					{	iii = jj;
#						olp = tlp;
#					}
#				}
#				ii = c(ii, iii);
#			}

			if(multi)
			{	r2 = summary(lm(yy~X[,ii]))$adj.r.squared;
			} else
			{	r2 = summary(lm(yy~apply(array(X[,ii],dim=c(nn,length(ii))),1,mean)))$adj.r.squared;
			}
			if(r2 < 0.01) next;
			dr2 = rep(r2, length(ii));
			if(length(ii)> 1)
			{	for(jj in 1:length(ii))
				{	if(multi)
					{	dr2[jj] = r2 - summary(lm(yy~X[,ii[-jj]]))$adj.r.squared;
					} else
					{	dr2[jj] = r2 - summary(lm(yy~apply(array(X[,ii[-jj]],dim=c(nn,length(ii)-1)),1,mean)))$adj.r.squared;
					}
				}
			}
			if(length(ii)>0) ii = rinv[j,1]+ii-1;
			ii = paste(ii,collapse=",");
			dr2 = paste(round(dr2*1000)/1000,collapse=",");
			rt = paste(paste(as.matrix(gene[gid[j],1:7]),collapse=" "),round(r2*1000)/1000,ii,dr2,collapse=" ");
			print(rt);

			write.table(rt, fout, quote=F,row.names=F,col.names=F,append=(j>1));
		
			rt = paste(as.matrix(gene[gid[j],1]),paste(rti$mcmc[,1]+rinv[j,1]-1,collapse=","),paste(round(rti$mcmc[,2]*1000)/10,collapse=","), collapse=" ");
			write.table(rt, paste(fout, ".mcmc", sep=""), quote=F,row.names=F,col.names=F,append=(j>1));
		}
	}
}

extractDistLoci<-function(fname)
{	x = read.table(fname, sep=" ");
	tss = x[,3];
	ort = x[,5];
	tss[ort<0] = x[ort<0,4];
	
	options(scipen=1000);
	chr = paste("chr",unique(as.matrix(x[,2])), sep="");
	library("data.table");
	region = NULL;
	for(i in chr)
	{	pos=as.matrix(fread(paste("/gpfs/group/yzz2/default/scratch/roadmap_analysis/tmp/test.",i,".state",sep=""))[,3]);
		tt = which(paste("chr", as.matrix(x[,2]), sep="")  == i);
		for(j in tt)
		{	l = as.integer(unlist(strsplit(as.matrix(x[j,9]),",")));
			dr2 = as.numeric(unlist(strsplit(as.matrix(x[j,10]),",")));
			md = (pos[l] - tss[j]) * ort[j];
			if(length(l)==1)
			{	#region = rbind(region, c(i, pos[l[ttt]]+100-800, pos[l[ttt]]+100+800, paste(j,"_",ttt,sep="")));
				region = rbind(region, c(j,x[j,8],l, md, dr2));
			} else
			{	#region = rbind(region, cbind(i, pos[l[ttt]]+100-800, pos[l[ttt]]+100+800, paste(j,"_",ttt,sep="")));
				region = rbind(region, cbind(j,x[j,8],l,md, dr2));
			}
			
			if(j%%100==0) print(j);
		}
	}
	write.table(region, paste(fname, ".bed", sep=""), quote=F, row.names=F, col.names=F);
#	mergeBED(paste(fname, ".bed", sep=""), paste(fname, ".bed", sep=""));
}

mergeBED<-function(fname, fout)
{	x=read.table(fname);
	chr = unique(as.matrix(x[,1]));
	nx = NULL;
	for(i in chr)
	{	t = which(x[,1]==i);
		nx = rbind(nx, x[t[order(x[t,2])],]);
	}

	chr = as.matrix(nx[,1]);
	pos = as.matrix(nx[,2:3]);
	id = as.matrix(nx[,4]);
	l = dim(nx)[1];
	ed = which(chr[-1] != chr[-l] | (chr[-1] == chr[-l] & pos[-1,1] > pos[-l,2]));
	st = c(1, ed + 1);
	ed = c(ed, l);
	
	nchr = chr[st];
	npos = cbind(pos[st,1], pos[ed,2]);
	nid = id[st];
	k = max(ed - st);
	if(k > 0)
	{	for(i in 1:k)
		{	st = st + 1;
			t = which(st <= ed);
			nid[t] = paste(nid[t], id[st[t]], sep=",");
		}
	}	

	write.table(cbind(nchr, npos, nid), fout, quote=F,row.names=F,col.names=F);
}

TFenrich<-function(fname, r2cut=0.5, drange = rbind(c(-1e6,1e6)), distance = 100, rmtss = FALSE, tssonly = 0, fnull = NULL, factorfile="wgEncodeRegTfbsClusteredV3.bed.gz", matchgene=FALSE)
{	library("data.table");
	gene = read.table(fname, sep=" ");
	if(length(fnull)==0)
	{	x = as.matrix(fread(paste(fname,".bed",sep="")));
	} else
	{	x = as.matrix(fread(fnull));
	}
	t = which(x[,2]>r2cut);
	m = as.matrix(x[,1])[t];
	d = as.matrix(x[,4])[t];
	x = as.matrix(x[,3])[t];

	chr = paste("chr",gene[1,2],sep="");
	g = fread(paste("/gpfs/group/yzz2/default/scratch/roadmap_analysis/test/test/test.",chr,".state",sep=""));
	pos = as.matrix(g[,3]);
	z = read.table(factorfile);
	tfz = sort(unique(z[,4]));
	z = z[z[,1]==chr,];

	xp=pos[x]/200;
	s=-sign(d);#-gene[m,5] * ( - sign(d));
	if(length(drange) == 0) 
	{	t = 1:length(x);
	} else
	{	t = NULL;
		drange = array(drange, dim=c(length(drange) / 2, 2));
		for(i in 1:dim(drange)[1])
		{	t = c(t, which(d>=drange[i,1] & d<=drange[i,2]));
		}
		t = unique(t);
	}

	tl = rep(0,max(pos)+100000);
	tss = gene[,3];
	tss[gene[,5]<0] = gene[gene[,5]<0,4];
	print(c(length(tss),length(tl)));
	if(rmtss | tssonly) 
	{ 	allgene=read.table("Ensembl_v65.Gencode_v10.ENSG.gene_info");
		allgene=allgene[allgene[,2]==substr(chr,4,100),];
		alltss = allgene[,3];
		alltss[which(allgene[,5]<0)] = allgene[which(allgene[,5]<0),4];
		if(length(tl)<max(alltss)+5000) tl = rep(0,max(alltss)+10000);
		for(i in alltss) tl[(i-5000):(i+5000)]=1; 
		print(mean(tl));
	}
	print(c(length(t),length(tl)));
	rr = NULL;
	for(i in 1:length(t))
	{	if(tl[xp[t[i]]*200+100]==1-tssonly) rr=c(rr,i);
	}
	if(length(rr)>0) t=t[-rr];
	print(c(length(rr),length(t)));

	mgene = as.matrix(gene[m,1]);
	ugene = unique(mgene);
	mgid = match(mgene, ugene);
	tgid = rep(-1,dim(z)[1]);
	if(matchgene)
	{	tgid = match(as.matrix(z[,5]), ugene);
	}

	t = t[order(xp[t])];
	p = round((z[,2]+z[,3])/2/200);
	tm = match(z[,4], tfz);
	cn = 0;
	st = 1;
	ed = 1;
	l = length(p);
	a = array(0,dim=c(length(tfz),distance * 2 + 1));
	apply(cbind(xp[t],s[t], mgid[t]),1,function(x){
		while(p[st] < x[1]-distance & st < l) { st<<-st+1; }
		while(p[ed] < x[1]+distance & ed < l) { ed<<-ed+1; }
		if(ed < l & ed > st) ed <<- ed - 1;
		k = (p[st:ed] - x[1])*x[2]+1+distance;
		o = tm[st:ed];
		tt = which(k > 0 & k <= 1+distance*2 & tl[p[st:ed]*200]==tssonly);
		if(matchgene) 
		{	ttt = which(tgid[st:ed] == x[3]);
			tt = intersect(ttt, tt);
		}
		k = k[tt];
		o = o[tt];
		if(length(k)>0)
		#{	apply(array(c(k,o),dim=c(length(k),2)),1,function(zz){a[zz[2],]<<-a[zz[2],]+1/(abs((1:(2*distance+1))-zz[1])+1);return(NULL)});
		{	apply(array(c(k,o),dim=c(length(k),2)),1,function(zz){a[zz[2],zz[1]]<<-a[zz[2],zz[1]]+1;return(NULL)});
		}
		if(cn%%100==0) print(c(cn,x[1],p[st],p[ed],l));
		cn <<- cn + 1;
		return(NULL);
	});

	n = tabulate(tm, nbins = length(tfz));
	rownames(a) = tfz;
	return(a);
}

stepMerge<-function(x, y, p, burnin = 50, mcmc = 50, step = 1, mmax = FALSE)
{	n = length(y);
	l = length(p);	
	if(length(x)/n==1)
		x = array(x,dim=c(n,length(x)/n));
	sel = c(rbinom(l, 1, p));
	k = sum(sel);
#	y = y - mean(y);
#	x = t(t(x) - apply(x,2,mean));
	xy = c(t(x)%*%y);
	dprior = log((p + 1e-5) / (1-p + 1e-5));

	lp = lpnull = lp0 = 0;
	X = rep(0, n);
	XY = 0;
	#v = max(sum(y^2)/(n-1),0.1);
	lpnull = getlp(X,y,XY,k,n);#(-sum(y^2)/v - n*log(v))/2;
	if(k == 0)
	{	lp0 = lp = lpnull;
	} else
	{	X = apply(array(x[,sel==1],dim=c(n,k)),1,sum);
		XY = sum(xy[sel==1]);
		#dy = y - XY / sum(X^2 + 1e-5) * X;
		#v = max(sum(dy^2)/(n-1),0.1);
		lp0 = getlp(X,y,XY,k,n);#(- n * (1 + log(v))) / 2;
		lp = lp0 + sum(log(p[sel==1]/(1-p[sel==1])));
	}
		
	s = rep(0, l);
	for(i in 1:(burnin + mcmc))
	{	#dy = y - XY / sum(X^2 + 1e-5) * X;
		#v = var(dy);
		#tlp = (-sum(dy^2)/v-n*log(v))/2;
		#if(sum(sel)>0)
		#{	tlp1 = logLik(lm(y~apply(array(x[,sel==1],dim=c(nn,sum(sel))),1,sum)-1));
	#		tlp1 = tlp1 + sum(log(p[sel==1]/(1-p[sel==1]))) + nn / 2 * log(2 * pi);
#			tlp = tlp + sum(log(p[sel==1]/(1-p[sel==1])));
#		} else
#		{	tlp1 = logLik(lm(y~1)) + nn / 2 * log(2 * pi); 
#		}
		#cat(i,"");#, tlp, tlp1));
		takemax = (i > burnin & mmax);
		o = sample(l);
		flag = !takemax;
		#flip
		u = runif(length(o));
		u = log(u / (1 - u));
		if(takemax) u = rep(0, length(o));
		#DLP = rep(-1000000,length(o));
		#for(j in o)
		apply(array(o,dim=c(length(o),1)),1,function(j)
		{	ts = 1 - 2 * sel[j];
			nX = X + ts * x[,j];
			nXY = XY + ts * xy[j];
			nk = k + ts;
			if(nk == 0)
			{	nlp = lpnull;
			} else
			{	#dy = y - nXY / sum(nX^2 + 1e-5) * nX;
				#v = max(sum(dy^2)/(n-1), 0.1);
				#nlp = (- n * (1 + log(v))) / 2;
				nlp = getlp(nX,y,nXY,nk,n);
			}
			dlp0 = nlp - lp0;
			dlp = dlp0 + ts * dprior[j];

			#DLP[j] <<- dlp;
			if(u[j] < dlp)
			{	X <<- nX;
				XY <<- nXY;
				sel[j] <<- 1 - sel[j];
				k <<- nk;
				lp0 <<- lp0 + dlp0;
				lp <<- lp + dlp;
				flag <<- TRUE;
			}
			return(NULL);
		});

		if(k > 0) #switch
		{	
			for(j in 1:(k * 2))
			{	if(k == 1)
				{	oi = which(sel==1)[1];
				} else
				{	oi = sample(which(sel==1),size=1);
				}
				cn = 0;
				while(cn < 20)
				{	cn = cn + 1;
					ni = oi + (rpois(1,5) + 1) * (2 * as.integer(runif(1)<0.5) - 1);
					if(ni > 0 & ni <= l)
						if(sel[ni] == 0) cn = 200;
				}
				if(cn < 200) next;
				nX = X - x[,oi] + x[,ni];
				nXY = XY - xy[oi] + xy[ni];
				#dy = y - nXY / sum(nX^2 + 1e-5) * nX;
				#v = max(sum(dy^2)/(n-1),0.1);
				dlp0 = getlp(nX,y,nXY,k,n) - lp0;#(- n * (1 + log(v))) / 2 - lp0;
				dlp = dlp0 + dprior[ni] - dprior[oi];
				u = runif(1);
				if(takemax) u = 0.5;
				if(log(u/(1-u)) < dlp)
				{	
					sel[oi] = 0; sel[ni] = 1;
					X = nX;
					XY = nXY;
					lp0 = lp0 + dlp0;
					lp = lp + dlp;
					flag = TRUE;
				}
			}
		}

		if(i > burnin & i%%step==0) s = s + sel;
		if(takemax) s = sel;
		if(takemax & !flag) break;
	}
	#cat("\n");
	return(s);
}

stepMulti<-function(x, y, p, burnin = 50, mcmc = 50, step = 1, mmax = FALSE)
{	n = length(y);
	l = length(p);	
	x = array(x,dim=c(n,length(x)/n));
	sel = c(rbinom(l, 1, p));
#	y = y - mean(y);
#	x = t(t(x) - apply(x,2,mean));
	dprior = log((p + 1e-5) / (1-p + 1e-5));

	lp0 = getlp_multi(x,y,which(sel==1),n);
	lp = lp0 + sum(log(p[sel==1]/(1-p[sel==1])));
		
	s = rep(0, l);
	for(i in 1:(burnin + mcmc))
	{	cat(i,"");
		takemax = (i > burnin & mmax);
		o = sample(l);
		flag = !takemax;
		#flip
		u = runif(length(o));
		if(takemax) u = rep(0.5, length(o));
		u = log(u / (1 - u));
		DLP = rep(-1000000,length(o));
		tt = which(sel == 1);
		apply(array(o,dim=c(length(o),1)),1,function(j)
		{	if(sel[j] == 1) {	ts = tt[-which(tt==j)];	}
			else { ts = c(tt, j);	}
			nlp = getlp_multi(x,y,ts,n);
			dlp0 = nlp - lp0;
			dlp = dlp0 + (1 - 2 * sel[j]) * dprior[j];

			DLP[j] <<- dlp;
			if(u[j] < dlp)
			{	tt <<- ts;
				sel[j] <<- 1 - sel[j];
				lp0 <<- lp0 + dlp0;
				lp <<- lp + dlp;
				flag <<- TRUE;
			}
			return(NULL);
		});

		tt = which(sel == 1);
		k = length(tt);
		if(k > 0) #switch
		{	
			for(j in 1:(k * 2))
			{	if(k == 1)
				{	oi = which(sel==1)[1];
				} else
				{	oi = sample(which(sel==1),size=1);
				}
				cn = 0;
				while(cn < 20)
				{	cn = cn + 1;
					ni = oi + (rpois(1,5) + 1) * (2 * as.integer(runif(1)<0.5) - 1);
					if(ni > 0 & ni <= l)
						if(sel[ni] == 0) cn = 200;
				}
				if(cn < 200) next;
				ts = tt[-which(tt==oi)];
				ts = c(tt, ni);
				dlp0 = getlp_multi(x,y,ts,n) - lp0;
				dlp = dlp0 + dprior[ni] - dprior[oi];
				u = runif(1);
				if(takemax) u = 0.5;
				if(log(u/(1-u)) < dlp)
				{	tt = ts;
					sel[oi] = 0; sel[ni] = 1;
					lp0 = lp0 + dlp0;
					lp = lp + dlp;
					flag = TRUE;
				}
			}
		}

		if(i > burnin & i%%step==0) s = s + sel;
		if(takemax) s = sel;
		if(takemax & !flag) break;
	}
	#cat("\n");
	return(s);
}

stepMerge_inc<-function(x, y, penalty, sel=NULL, sz = 0)
{	n = length(y);
	x = array(x, dim=c(n, length(x)/n));
	l = dim(x)[2];	
	#y = y - mean(y);
	#x = t(t(x) - apply(x, 2, mean));
	xy = c(t(x)%*%y);

	lp = lpnull = 0;
	X = rep(0, n);
	XY = 0;
	#v = var(y);
	sk = length(sel);
	lpnull = getlp(X,y,XY,sk,n);#(-sum(y^2)/v - n*log(v))/2;
	if(sk == 0)
	{	lp = lpnull;
	} else
	{	X = apply(array(x[,sel],dim=c(n,sk)),1,sum);
		XY = sum(xy[sel]);
		#dy = y - XY / sum(X^2 + 1e-5) * X;
		#v = var(dy);
		lp = getlp(X,y,XY,sk,n);#(-sum(dy^2)/v - n * log(v)) / 2;
	}
	
	for(i in 1:l)
	{	j = 1;
		nlp = apply(x,2,function(z)
			{	#dy = y-(XY+xy[j])/sum((X+z)^2+1e-5)*(X+z);
				#v = var(dy);
				tlp = getlp(X+z,y,XY+xy[j],sk+1,n);
				j <<- j+1;
				#return((-sum(dy^2)/v-n*log(v))/2);
				return(tlp);
			});
		if(length(sel) > 0) nlp[sel] = lp-1000000;
		if(sz > 0 | max(nlp) > lp + penalty)
		{	k = which(nlp == max(nlp))[1];
			sel = c(sel, k);
			sk = sk + 1;
			X = X + x[,k];
			XY = XY + xy[k];
			lp = max(nlp);
			#cat(k,"");
		} else
		{	break;
		}
		if(sz > 0 & i >= sz) break;
	}
	#cat("\n");
	return(sel);
}

stepMerge_dec<-function(x, y, penalty, sel=NULL)
{	n = length(y);
	x = array(x, dim=c(n, length(x)/n));
	l = dim(x)[2];	
#	y = y - mean(y);
#	x = t(t(x) - apply(x,2,mean));
	xy = c(t(x)%*%y);

	lp = lpnull = 0;
	X = rep(0, n);
	XY = 0;
	sel = 1:l;
	sk = length(sel);
#	v = var(y);
	lpnull = getlp(X,y,XY,sk,n);#(-sum(y^2)/v - n*log(v))/2;
	if(sk == 0)
	{	lp = lpnull;
	} else
	{	X = apply(array(x[,sel],dim=c(n,sk)),1,sum);
		XY = sum(xy[sel]);
		#dy = y - XY / sum(X^2 + 1e-5) * X;
		#v = var(dy);
		lp = getlp(X,y,XY,sk,n);#(-sum(dy^2)/v - n * log(v)) / 2;
	}
	
	for(i in 1:length(sel))
	{	j = 1;
		nlp = apply(array(x[,sel],dim=c(n,length(sel))),2,function(z)
			{	#dy = y-(XY-xy[sel[j]])/sum((X-z)^2+1e-5)*(X-z);
				#v = var(dy);
				tlp = getlp(X-z,y,XY-xy[sel[j]],sk-1,n);
				j <<- j+1;
				#return((-sum(dy^2)/v-n*log(v))/2);
				return(tlp);
			});
		if(max(nlp) > lp - penalty)
		{	k = which(nlp == max(nlp))[1];
			X = X - x[,sel[k]];
			XY = XY - xy[sel[k]];
			lp = max(nlp);
			sel = sel[-k];
			sk = sk - 1;
			#cat(k,"");
		} else
		{	
			break;
		}
	}
	#cat("\n");
	return(sel);
}

sampleNull<-function(fcase, chr, path, statepref, maxdist = 200000, fpref = "stepnull")
{	case = read.table(fcase,sep=" ");
	casebed = read.table(paste(fcase,".bed",sep=""));

	gene=read.table("Ensembl_v65.Gencode_v10.ENSG.gene_info");
	Y=array(0,dim=c(dim(gene)[1],56));
	x=read.table("57epigenomes.RPKM.pc",header=T);
	cells=colnames(x)[-(1:2)];
	m=match(gene[,1],x[,1]);
	t=which(is.na(m)==F);
	Y[t,]=as.matrix(x[m[t],-(1:2)]);
	x=read.table("57epigenomes.RPKM.nc",header=T);
	m=match(gene[,1],x[,1]);
	t=which(is.na(m)==F);
	Y[t,]=as.matrix(x[m[t],-(1:2)]);
	x=read.table("57epigenomes.RPKM.rb",header=T);
	m=match(gene[,1],x[,1]);
	t=which(is.na(m)==F);
	Y[t,]=as.matrix(x[m[t],-(1:2)]);

	Y=log2(Y+0.01);

	t=which(apply(Y,1,sd) > 0);
	Y=Y[t,];
	gene=gene[t,];	
	genechr = paste("chr",gene[,2],sep="");
	tss = gene[,3];
	tss[gene[,5]<0] = gene[gene[,5]<0,4];
	ort = gene[,5];

	library("data.table");
	library("MASS");
	cellinfo = read.table("/gpfs/group/yzz2/default/RoadmapEpigenomics/extra/epigenomeID.txt",header=T,comment="!",sep="\t");
	cellinfo = cellinfo[order(as.integer(substr(cellinfo[,1],2,1000))),];
	cellmap = match(cells, cellinfo[,1]);
	nn = dim(Y)[2];
	for(i in chr)
	{	print(i);
		g = fread(paste(path,statepref,i,".state",sep=""));
		a=system(paste("wc -l ",path, statepref,i,".para",sep=""),intern=T);
		gn=as.integer(substr(a,1,regexpr(" ",a)[1]))-1;
		pos = as.matrix(g[,3]);
		state = as.matrix(g[,5:131])[,cellmap];
		t = which(gene[,2]==substr(i,4,1000));

		print("state file loaded");
	
		gid = tinv = rinv = NULL;
		for(j in t)
		{	tt = which(pos >= tss[j] - 1000 & pos < tss[j] + 1000);
			if(length(tt) == 0) next;
			rr = which(pos >= tss[j] - maxdist & pos < tss[j] + maxdist);
			if(length(rr) < 2) next;
			gid = c(gid, j);
			tinv = rbind(tinv, range(tt));
			rinv = rbind(rinv, range(rr));	
		}

		if(length(gid) < 10) next;
		#e = c(as.matrix(lm(c(Y[gid,])~as.factor(state[as.integer((tinv[,1]+tinv[,2])/2),])-1)$coef));
		#e[is.na(e)==T]=0;
		ee = NULL;
		str = paste(path,statepref,"coef",sep="");
		if(file.exists(str) == F)
		{	ee = getStateCoef(paste(path,statepref,sep=""));
		} else
		{	ee = read.table(str);
		}
		
		lab = rownames(ee);
		l = apply(t(t(lab)),1,function(x){tail(unlist(gregexpr("_",x)),n=1)});
		lab = substr(lab,1,l-1);
		ee = t(array(ee[,(dim(ee)[2]+1)/2],dim=c(gn,dim(ee)[1]/gn)));
		rownames(ee) = unique(lab);
		colnames(ee) = 1:gn-1;
		t=which(apply(ee,1,sd) == 0);
		if(length(t)>0) ee[t,] = ee[dim(ee)[1],];

		fout = paste("/gpfs/scratch/yzz2/rna/", fpref, ".", statepref,i,".bed",sep="");
		fout0 = paste("/gpfs/scratch/yzz2/rna/", fpref, ".", statepref,i,".bed0",sep="");
		sn = rep(0,length(gid));
		rr = rrb = NULL;

		kk = 1;
		tindex = NULL;
		for(j in 1:max(casebed[,1]))
		{	tt=which(casebed[,1]==j);
			while(as.matrix(gene[gid[kk],1])!=as.matrix(case[j,1])) kk=kk+1;
			tindex = c(tindex, rep(kk, length(tt)));
			kk = kk + 1;
		}
		uit = unique(tindex);
		for(j in uit)
		{	
			if(length(ee)>0) e = ee[which(gene[gid[j],6]==rownames(ee)),];	
			X = t(array(e[c(state[rinv[j,1]:rinv[j,2],])+1], dim=c(rinv[j,2]-rinv[j,1]+1,nn)));
			#X = smoothMatrix(X);
			yy = Y[gid[j],];
			r = cor(yy, X+rnorm(length(X),0,0.0001));
			r[is.na(r)==T] = 0;
			#a = 1:dim(X)[2];
			a = which(r > minr);
			tt=which(uit==j);
			tt=which(casebed[,1]==tt);
			rr = c(rr, r[casebed[tt,3]-rinv[j,1]+1]);
			rrb = c(rrb, r[a]);
			sn[j] = length(tt);
			print(c(j, length(a),length(tt),sn[j]));
		}
		tn = (density(rr,from=minr,to=1,n=100)$y+0.01)/(density(rrb,from=minr,to=1,n=100)$y+0.01);
		#tn = (density(rr,from=-1,to=1,n=100)$y+0.01)/(density(rrb,from=-1,to=1,n=100)$y+0.01);
		for(j in uit)
		{	X = t(array(e[c(state[rinv[j,1]:rinv[j,2],])+1], dim=c(rinv[j,2]-rinv[j,1]+1,nn)));
			#X = smoothMatrix(X);
			yy = Y[gid[j],];
			r = cor(yy, X+rnorm(length(X),0,0.0001));
			r[is.na(r)==T] = 0;
			#a = 1:dim(X)[2];
			#w = round((r[a]--1)/2*100)+1;
			a = which(r > minr);
			w = round((r[a]-minr)/(1-minr)*100)+1;
			w[w>100]=100;
			w = tn[w];
			if(length(a)<2) 
			{	t0 = a;
				t00 = a;
			} else
			{	t0 = sort(sample(a,size=sn[j],prob=w,replace=T));
				t00 = sort(sample(a,size=sn[j],replace=F));
			}
			tt=which(uit==j);
			tt=which(casebed[,1]==tt);
			md = (pos[t0 + rinv[j,1] - 1] - tss[gid[j]]) * ort[gid[j]];
			md0 = (pos[t00 + rinv[j,1] - 1] - tss[gid[j]]) * ort[gid[j]];
			if(length(tt)==1)
			{	rt = rt0 = casebed[tt,];
				rt[3] = t0 + rinv[j,1] - 1;
				rt0[3] = t00 + rinv[j,1] - 1;
				rt[4] = md;
				rt0[4] = md0;
				rt = c(rt, round(r[casebed[tt,3]-rinv[j,1]+1]*1000)/1000,round(r[t0]*1000)/1000);
				rt0 = c(rt0, round(r[casebed[tt,3]-rinv[j,1]+1]*1000)/1000,round(r[t00]*1000)/1000);
			} else
			{	rt = rt0 = casebed[tt,];
				rt[,3] = t0 + rinv[j,1] - 1;
				rt0[,3] = t00 + rinv[j,1] - 1;
				rt[,4] = md;
				rt0[,4] = md0;
				rt = cbind(rt, round(r[casebed[tt,3]-rinv[j,1]+1]*1000)/1000,round(r[t0]*1000)/1000);
				rt0 = cbind(rt0, round(r[casebed[tt,3]-rinv[j,1]+1]*1000)/1000,round(r[t00]*1000)/1000);
			}
			cat(j,"");
			write.table(rt, fout, quote=F,row.names=F,col.names=F,append=(j>uit[1]));
			write.table(rt0, fout0, quote=F,row.names=F,col.names=F,append=(j>uit[1]));
		}
	}
}

runStepMerge<-function(x, y, p = NULL, keep=50, B = 20, burnin = 100, mcmc = 400, direction="inc", penalty = 1)
{	l = length(x) / length(y);
	x = array(x,dim=c(length(y),l));
	ss=rep(0,l);
	if(length(p)==0) p = rep(0.1,l);#min(0.5, 20 / l), l);
	for(i in 1:B)
	{	tss=stepMerge(x,y,p,burnin=burnin/B,mcmc=mcmc/B);
		if(i==1){ss=tss;}else{ss=ss+tss}
	}
	tt=order(ss+c(cor(x+rnorm(length(x),0,0.000001),y))/100)[l+-(min(keep,l)-1):0];
	tt = tt[which(ss[tt]>1)];
	rt=NULL;
	rt$mcmc = rt$ss = NULL;
	if(length(tt)==0) return(rt);

	rt$mcmc=cbind(tt,ss[tt]/mcmc);
	ss=rep(0,l);
	tp = p[tt];
	tp = 1-(1-tp)^(l/length(tt));
	for(i in 1:B)
	{	ss[tt]=ss[tt]+stepMerge(x[,tt],y,tp,burnin=burnin/B,mcmc=mcmc/B);
	}
	if(max(ss)/mcmc< 0.05) ss[which(ss==max(ss))]=mcmc;
	tt = which(ss/mcmc >= 0.05);
print(length(tt));

	if(length(tt)==0)
	{	rt$ss = NULL;
		return(rt);
	}
	if(direction == "inc")
	{	ss = sort(tt[stepMerge_inc(x[,tt],y,sel=which(p[tt]>0.99),penalty=penalty)]);
	} else
	{	ss = sort(tt[stepMerge_dec(x[,tt],y,penalty=penalty)]);#, sel=which(p[tt]>0.99))]);
	}

	rt$ss = ss;

	return(rt);
}

DEpredict<-function(chr, path, statepref, maxdist = 200000, fpref = "step", direction="inc", gidlist = NULL)
{	gene=read.table("Ensembl_v65.Gencode_v10.ENSG.gene_info");
	Y=array(0,dim=c(dim(gene)[1],56));
	x=read.table("57epigenomes.RPKM.pc",header=T);
	cells=colnames(x)[-(1:2)];
	m=match(gene[,1],x[,1]);
	t=which(is.na(m)==F);
	Y[t,]=as.matrix(x[m[t],-(1:2)]);
	x=read.table("57epigenomes.RPKM.nc",header=T);
	m=match(gene[,1],x[,1]);
	t=which(is.na(m)==F);
	Y[t,]=as.matrix(x[m[t],-(1:2)]);
	x=read.table("57epigenomes.RPKM.rb",header=T);
	m=match(gene[,1],x[,1]);
	t=which(is.na(m)==F);
	Y[t,]=as.matrix(x[m[t],-(1:2)]);

	Y=log2(Y+0.01);

	t=which(apply(Y,1,sd) > 0);
	Y=Y[t,];
	gene=gene[t,];	
	genechr = paste("chr",gene[,2],sep="");
	tss = gene[,3];
	tss[gene[,5]<0] = gene[gene[,5]<0,4];

	library("MASS");
	cellinfo = read.table("/gpfs/group/yzz2/default/RoadmapEpigenomics/extra/epigenomeID.txt",header=T,comment="!",sep="\t");
	cellinfo = cellinfo[order(as.integer(substr(cellinfo[,1],2,1000))),];
	cellmap = match(cells, cellinfo[,1]);
	nn = dim(Y)[2];
	
	library("data.table");
	i = chr;
	print(i);
	g = fread(paste(path,statepref,i,".state",sep=""));
	a=system(paste("wc -l ",path, statepref,i,".para",sep=""),intern=T);
	gn=as.integer(substr(a,1,regexpr(" ",a)[1]))-1;
	pos = as.matrix(g[,3]);
	state = as.matrix(g[,5:131])[,cellmap];
	t = which(gene[,2]==substr(i,4,1000));
	
	gid = tinv = rinv = NULL;	
	for(j in t)
	{	tt = which(pos >= tss[j] - 1000 & pos < tss[j] + 1000);
		if(length(tt) == 0) next;
		rr = which(pos >= tss[j] - maxdist & pos < tss[j] + maxdist);
		if(length(rr) < 2) next;
		gid = c(gid, j);
		tinv = rbind(tinv, range(tt));
		rinv = rbind(rinv, range(rr));	
	}

	#e = c(as.matrix(lm(c(Y[gid,])~as.factor(state[as.integer((tinv[,1]+tinv[,2])/2),])-1)$coef));
	#e[is.na(e)==T]=0;
	ee = NULL;
	str = paste(path,statepref,"coef",sep="");
	if(file.exists(str) == F)
	{	ee = getStateCoef(paste(path,statepref,sep=""));
	} else
	{	ee = read.table(str);
	}
	
	lab = rownames(ee);
	l = apply(t(t(lab)),1,function(x){tail(unlist(gregexpr("_",x)),n=1)});
	lab = substr(lab,1,l-1);
	ee = t(array(ee[,(dim(ee)[2]+1)/2],dim=c(gn,dim(ee)[1]/gn)));
	rownames(ee) = unique(lab);
	colnames(ee) = 1:gn-1;
	t=which(apply(ee,1,sd) == 0);
	if(length(t)>0) ee[t,] = ee[dim(ee)[1],];

	fout = paste("/gpfs/scratch/yzz2/rna/", fpref, ".", statepref,i,".predict",sep="");
	fout1 = paste("/gpfs/scratch/yzz2/rna/", fpref, ".", statepref,i,".predict_crelist",sep="");
	if(length(gidlist) == 0) gidlist = 1:length(gid);
	gidlist = gidlist[which(gidlist>0 & gidlist<=length(gid))];
	for(j in gidlist)
	{	
		if(length(ee)>0) e = ee[which(gene[gid[j],6]==rownames(ee)),];	
		X = t(array(e[c(state[rinv[j,1]:rinv[j,2],])+1], dim=c(rinv[j,2]-rinv[j,1]+1,nn)));
		#X = smoothMatrix(X);
		yy = Y[gid[j],];
		
		tl = NULL;
		ff = NULL;
		iii = NULL;
		for(k in 1:nn)
		{	ff[k] = mean(yy[-k]);
			r = cor(yy[-k], X[-k,]+rnorm(length(X[-k,]),0,0.0001));
			r[is.na(r)==T] = 0;
			#a = 1:dim(X)[2];
			a = which(r > minr);
			if(length(a) < 1) next;
			#rti = runStepMerge(X[-k,a], yy[-k], mcmc=100, direction = direction);#log(nn-1)/4);
			rti = runStepMerge1(X[-k,a], yy[-k], mcmc=100, direction = direction,penalty=0.0);#log(nn-1)/4);
			ii = a[rti$ss];
			#ii=a[rti$mcmc[order(rti$mcmc[,2])[dim(rti$mcmc)[1]+-ceiling(sum(rti$mcmc[,2])):0],1]]
			if(length(ii) < 1) next;
			tl = c(tl,length(ii));

			tx = apply(array(X[,ii],dim=c(nn,length(ii))),1,mean);
			rr = getbeta(tx[-k], yy[-k]);
			ff[k] = rr[1] + tx[k] * rr[2];
			message(k," (",round(yy[k]*100)/100,", ",round(ff[k]*100)/100,") ",paste(sort(ii),collapse=" "));
			ii = rinv[j,1]+ii-1;
			iii = c(iii, ii);
		}
		tll = table(sort(iii));
		tll = paste(names(tll),as.integer(tll),sep=":");
		boxplot(tl);print(c(j,mean(tl),cor(yy,ff)));
		flag = file.exists(fout);
		write.table(paste(c(as.matrix(gene[gid[j],1]),round(ff*100)/100),collapse=" "),fout, quote=F,row.names=F,col.names=F,append=flag); 
		flag = file.exists(fout1);
		write.table(paste(c(as.matrix(gene[gid[j],1]),round(cor(yy,ff)*1000)/1000,mean(tl),c(paste(tll,collapse=","))),collapse=" "),fout1, quote=F,row.names=F,col.names=F,append=flag); 
	}
}

showFactorExp<-function(mcmc, select, state, cellmap, lab, gene, gid,e=NULL,showe=F)
{
	tg=as.numeric(gene[match(mcmc[gid,1],as.matrix(gene[,1])),-(1:2)]);
	tg=log2(tg+0.1)*3;
	mg=mean(tg);
	sg=sd(tg);
tg=tg-mean(tg);
	l=as.integer(unlist(strsplit(as.matrix(mcmc[gid,2]),",")));
	p=as.numeric(unlist(strsplit(as.matrix(mcmc[gid,3]),",")));
	if(max(p) > 100) p = p / 50;
	o=order(l);
	l=l[o];
	p=p[o];
	d=20;
	mycol=c(as.matrix(lab[state[l,cellmap]+1,6]));
	defaultpal=NULL;
	if(length(e)>0 & showe)
	{	colors=0:100/100*(range(e)[2]-range(e)[1])+range(e)[1];
        	my_palette=colorRampPalette(c("cyan","blue","black","red","yellow"))(n=100);
		defaultpal = palette(my_palette);	
		mycol=ceiling((e[state[l,cellmap]+1]-range(e)[1])*100/(range(e)[2]-range(e)[1]));
		mycol[mycol==0]=1;
	}

	rr=max(tg)-min(0,min(tg));
	plot(-10,-10,xlim=c(0,length(l)+rr),ylim=c(0,56+d));
	rect(rep(1:length(l)-1,56),rep(0:55,each=length(l))+d,rep(1:length(l),56)-0.05,rep(1:56,each=length(l))-0.05+d,col=mycol,border="black");
	points(tg+max(0,-min(tg))+1+length(l),1:56-0.5+d,pch=5,cex=0.6);
	lines(rep(length(l)+1+max(0,-min(tg)),2),c(d,d+56));
	t=match(as.integer(unlist(strsplit(as.matrix(select[gid,9]),","))),l);
	s=rep(1,length(l));s[t]=2;
	lines(1:length(l)-0.5,p/100*d,type="h",col=c("blue","red")[s],lwd=2)

	rt=NULL;rt$x=rt$y=NULL;
	r2 = r2s = 0;
	if(length(e)>0)
	{	X = array(e[state[l,cellmap]+1],dim=c(length(l),56));
		tx = apply(array(X,dim=c(length(l),56)),2,mean);
		txs = apply(array(X[t,],dim=c(length(t),56)),2,mean);
		r2=summary(lm(tg~tx))$adj.r.squared;
		r2s=summary(lm(tg~txs))$adj.r.squared;
		rt$x=txs;rt$y=tg;
		rt$X=t(X[t,]);
		rect(length(l),-0.5,length(l)+rr+1,d);
		points((txs-min(txs))/(max(txs)-min(txs))*rr+length(l)+0.5,0.5+(tg-min(tg))*(d-0.5)/(max(tg)-min(tg)),pch=20,cex=0.5,col="blue");
		lines(c(length(l),length(l)+rr+1),c(-0.5,d));
		text(length(l)+rr-2,1,round(r2s*1000)/1000);
	}
	print(c(mean(tg),sd(tg),y[gid,8],r2, r2s));
	if(length(defaultpal)>0) palette(defaultpal);
	return(rt);
}

getlp<-function(x, y, XY, p, n, da=1, db=1)
{	f = NULL;
	if(p>0)
	{	x=x/p;
		XY=XY/p;
	}
	#x2 = sum(x^2);
#	y = y - x;
#	XY = XY - x2;
	nn = n + da;
	#X = sum(x);
	#X2 = x2 + db;
	#Y = sum(y);

	ty = sum(y) - y;
	tx2 = sum(x^2) + db - x^2;
	tx = sum(x) - x;
	txy = XY - x*y;# + x^2;
	b = cbind(ty * tx2 - tx * txy, - ty * tx + (nn - 1) * txy) / (tx2 * (nn - 1) - tx^2);
	f = c(apply(b*cbind(1,x),1,sum));#b[,1] + b[,2] * x;

#	for(i in 1:n)
#	{	ty = Y - y[i];
#		tx2 = X2 - x[i]^2;
#		tx = X - x[i];	
#		txy = XY - x[i]*y[i];# + x[i]^2;
#		b = c(ty * tx2 - tx * txy, - ty * tx + (nn - 1) * txy) / (tx2 * (nn - 1) - tx^2);
#		f[i] = b[1] + b[2] * x[i];
		#print(c(b,f[i]));
#	}
#print(cor(y,f));
	dy = y - f;
	ss = sum(dy^2);
	v = max(ss / (n - 1), 0.1);
	lp = -(ss / v + n * log(v)) / 2;
	return(lp);
}

getlp_core<-function(x, y, XY, p, n, da=1, db=1)
{	if(p>0)
	{	x = x / p;
		XY = XY / p;
	}
	x2 = sum(x^2);
	y = y - x;
	XY = XY - x2;

	nn = n + da;
	X = sum(x);
	X2 = x2 + db;
	Y = sum(y);
		
	beta = c(Y * X2 - X * XY, - Y * X + nn * XY) / (X2 * nn - X^2);
	dy = y - beta[1] - beta[2] * x;
	ss = sum(dy^2);
	v = max(ss / (n - 1), 0.1);
	lp = -(ss / v + n * log(v)) / 2;
	return(lp);
}

getlp_multi<-function(x, y, sel, n, da=1, db=1)
{	p = length(sel);
	tx = rep(1, n);
	if(p>0)
	{	y = y - apply(array(x[,sel],dim=c(n,p)), 1, mean);
		tx = cbind(tx, x[,sel] / p);
	}
	beta = solve(t(tx)%*%tx+diag(c(da,rep(db,p))))%*%t(tx)%*%y;
	f = c(tx%*%beta);
	dy = y - f;
	ss = sum(dy^2);
	v = max(ss / (n - 1), 0.1);
	lp = -(ss / v + n * log(v)) / 2;
	return(lp);
}

getbeta<-function(x,y,da=1,db=1)
{	x2 = sum(x^2);
	y = y - x;
	XY = sum(x * y);

	nn = length(x) + da;
	X = sum(x);
	X2 = x2 + db;
	Y = sum(y);
		
	beta = c(Y * X2 - X * XY, - Y * X + nn * XY) / (X2 * nn - X^2);
	beta[2] = beta[2] + 1;
	return(beta);
}

eQTLs<-function(fname, r2cut = 0.5, drange = rbind(c(-1e6,1e6)), distance = 250, rmtss = FALSE, tssonly = 0, fnull = NULL)
{	library("data.table");
	gene = as.matrix(fread("/gpfs/group/yzz2/default/scratch/roadmap_analysis/gtex/gene_id.txt",header=F));
	l=as.integer(regexpr("\\.",gene[,1]));
	gene = substr(gene[,1],1,l-1);
	
	snp = as.matrix(fread("/gpfs/group/yzz2/default/scratch/roadmap_analysis/gtex/var_id.txt",header=F));
	snpchr = as.matrix(snp[,1]);
	snppos = round(as.integer(as.matrix(snp[,2]))/200);	

	pair = as.matrix(fread("/gpfs/group/yzz2/default/scratch/roadmap_analysis/gtex/allsignif_snpgene_pairs.txt"));
	pair[,4] = -log10(pair[,4]);

	x = read.table(fname, sep=" ");
	if(length(fnull)==0)
	{	xb = as.matrix(fread(paste(fname,".bed",sep="")));
	} else
	{	xb = as.matrix(fread(fnull));
	}
	t = which(xb[,2]>r2cut);
	m = as.matrix(xb[,1])[t];
	d = as.matrix(xb[,4])[t];
	xb = as.matrix(xb[,3])[t];
	chr = as.character(as.matrix(x[1,2]));

	gid = match(as.matrix(x[m,1]),gene);
	t = which(snpchr[pair[,2]] == chr & is.na(match(pair[,3],gid[is.na(gid)==F]))==F);
	pair = pair[t,];
	pair[,4] = pair[,4]-5;
	pair = pair[pair[,4]>0,];

	g = fread(paste("/gpfs/group/yzz2/default/scratch/roadmap_analysis/test/test/test.chr",chr,".state",sep=""));
	pos = as.matrix(g[,3]);
	xp=pos[xb]/200;
	s=-sign(d);#x[m,5];# * ( - sign(d));
	if(length(drange) == 0) 
	{	t = 1:length(xb);
	} else
	{	t = NULL;
		drange = array(drange, dim=c(length(drange) / 2, 2));
		for(i in 1:dim(drange)[1])
		{	t = c(t, which(d>=drange[i,1] & d<=drange[i,2]));
		}
		t = unique(t);
	}

	tl = rep(0,max(pos)+100000);
	tss = x[,3];
	tss[x[,5]<0] = x[x[,5]<0,4];
	tss = as.integer(tss[m]/200);
#	xp=tss+round(runif(length(xp),-999,999));
#	xp[xp<1]=1;xp[xp>dim(g)[1]]=dim(g)[1];
	if(rmtss | tssonly) 
	{ 	allgene=read.table("Ensembl_v65.Gencode_v10.ENSG.gene_info");
		allgene=allgene[allgene[,2]==chr,];
		alltss = allgene[,3];
		alltss[which(allgene[,5]<0)] = allgene[which(allgene[,5]<0),4];
		if(length(tl)<max(alltss)+5000) tl = rep(0,max(alltss)+10000);
		for(i in alltss) tl[(i-5000):(i+5000)]=1; 
		print(mean(tl));
	}
	rr = which(is.na(gid[t])==T);
	for(i in 1:length(t))
	{	if(tl[xp[t[i]]*200+100]==1-tssonly) rr=c(rr,i);
	}
	if(length(rr)>0) t=t[-rr];
	t = t[order(xp[t])];

	cn = 0;
	a = b = array(0, dim=c(44,distance*2+1));
	apply(cbind(xp[t],s[t],gid[t],tss[t]),1,function(z){
		tt = which(pair[,3]==z[3]);
		if(length(tt)>0)
		{	l = snppos[pair[tt,2]];
			p = pair[tt,4];
			o = pair[tt,1];
			k = (l - z[1])*z[2]+1+distance;
			ttt = which(k > 0 & k <= 1+distance*2 & tl[l*200]==tssonly);
			#ttt = which(k > -distance & k <= 1+distance*3 & tl[l*200]==tssonly);
			if(length(ttt)>0)
			{	k = k[ttt];
				o = o[ttt];
				p = p[ttt];
				aa = bb = array(0, dim=c(44,1 + 2 * distance));
				aa[(k-1)*44+o] = p;
				bb[(k-1)*44+o] = 1;
				#apply(array(c(k,o,p),dim=c(length(k),3)),1,function(zz){a[zz[2],zz[1]]<<-a[zz[2],zz[1]]+zz[3];b[zz[2],zz[1]]<<- b[zz[2],zz[1]]+1;return(NULL)});
			#	apply(array(c(k,o,p),dim=c(length(k),3)),1,function(zz){dd=1/(abs(1:(2*distance+1)-zz[1])+1)^1.5; tdd=dd*zz[3] - aa[zz[2],];tdd[tdd<0]=0;aa[zz[2],]<<-aa[zz[2],]+tdd;tdd=dd-bb[zz[2],];tdd[tdd<0]=0;bb[zz[2],]<<-bb[zz[2],]+tdd;return(NULL)});
				a <<- a + aa;
				b <<- b + bb;
			}
			if(cn%%100==0) print(cn);
			cn <<- cn + 1;
		}
		return(NULL);
	});
#	a = a / length(t);
#	b = b / length(t);
	rt=NULL;
	rt$a=a;
	rt$b=b;
	return(rt);
}

showEQTL<-function(gid, select, statepos, gene, snppos, pair,snull=NULL)
{	tpos = statepos[as.integer(unlist(strsplit(as.matrix(select[gid,9]),",")))];
	tgene = as.matrix(select[gid,1]);
	m = match(tgene, gene);
	if(is.na(m)==T) return(NULL);
	t=which(pair[,3]==m);
	if(length(t)==0) return(NULL);
	epos = snppos[pair[t,2]]*200;
	eqval = pair[t,4];
	tss = select[gid,3];
	if(select[gid,5]<0) tss=select[gid,4];
	ort=select[gid,5];
	
	plot((epos-tss)*ort, eqval,pch=rep(20,length(epos)),type="h",ylim=c(0,max(eqval)),xlim=c(-1,1)*max(2e5,max(abs(epos-tss))),col=rep("black",length(epos)));
	points((tpos-tss)*ort,rep(0,length(tpos)),pch=17,col="red");
	if(length(snull)>0)
	{	tpos0=statepos[snull[snull[,1]==gid,3]];
		points((tpos0-tss)*ort,rep(0,length(tpos0)),pch=2,col="green");
	}
	rect(-200000,-0.2,200000,-0.1,col="blue",border=NA);
	rect(-5000,-0.2,5000,-0.1,col="cyan",border=NA);
}
 
runStepMerge1<-function(x, y, penalty = 1, B = 10, burnin = 100, mcmc = 400, direction = "inc")
{	l = length(x) / length(y);
	x = array(x,dim=c(length(y),l));
	sz = min(l, 40);
	nn = ceiling(l / sz);
	sz = ceiling(l / nn);
	ss = NULL;
	for(i in 1:nn)
	{	a = (i-1)*sz + 1:sz;
		a = a[a<=l];
tss=stepMerge(x[,a],y,rep(0.1,length(a)),burnin=1,mcmc=5);
a=a[which(tss>1)];
if(length(a)==0) next;
		#if(direction=="inc")
		#{	tss=stepMerge_inc(x[,a],y,penalty);
		#} else
		#{	tss=stepMerge_dec(x[,a],y,penalty);
		#}
		#tss=a[tss]-(i-1)*sz;
tss=a-(i-1)*sz;
		ss = c(ss, (i-1)*sz + tss);
	}
ss = unique(sort(ss));
#	ss = unique(c(ss, which(stepMerge(x,y,rep(0.1/nn,dim(x)[2]))>mcmc*0.1/nn*2)));

	rt=NULL;
	rt$mcmc = rt$ss = NULL;
	if(length(ss)==0) return(rt);

#print("b");
#	tss = stepMerge(x[,ss],y,rep(0.1,length(ss)),burnin = burnin, mcmc = mcmc);
#	tt = order(tss+cor(x[,ss]+rnorm(length(x[,ss]),0,0.0000001),y)/100)[length(ss)+-(min(length(ss),50)-1):0];
#	if(length(tt)==0) return(rt);
#	rt$mcmc = cbind(ss[tt], tss[tt]/mcmc);
#	ss = ss[tt];

	if(direction == "inc")
	{	ss = sort(ss[stepMerge_inc(x[,ss],y,penalty=penalty)]);
	} else
	{	ss = sort(ss[stepMerge_dec(x[,ss],y,penalty=penalty)]);#, sel=which(p[tt]>0.99))]);
	}
	
	rt$ss = ss;	
	return(rt);
}

check<-function(chr, path, statepref, maxdist=200000, fpref="step", pchic=FALSE)
{	gene=read.table("Ensembl_v65.Gencode_v10.ENSG.gene_info");
	Y=array(0,dim=c(dim(gene)[1],56));
	x=read.table("57epigenomes.RPKM.pc",header=T);
	cells=colnames(x)[-(1:2)];
	m=match(gene[,1],x[,1]);
	t=which(is.na(m)==F);
	Y[t,]=as.matrix(x[m[t],-(1:2)]);
	x=read.table("57epigenomes.RPKM.nc",header=T);
	m=match(gene[,1],x[,1]);
	t=which(is.na(m)==F);
	Y[t,]=as.matrix(x[m[t],-(1:2)]);
	x=read.table("57epigenomes.RPKM.rb",header=T);
	m=match(gene[,1],x[,1]);
	t=which(is.na(m)==F);
	Y[t,]=as.matrix(x[m[t],-(1:2)]);

	Y=log2(Y+0.01);

	t=which(apply(Y,1,sd) > 0);
	Y=Y[t,];
	gene=gene[t,];	
	genechr = paste("chr",gene[,2],sep="");
	tss = gene[,3];
	tss[gene[,5]<0] = gene[gene[,5]<0,4];

	library("MASS");
	cellinfo = read.table("/gpfs/group/yzz2/default/RoadmapEpigenomics/extra/epigenomeID.txt",header=T,comment="!",sep="\t");
	cellinfo = cellinfo[order(as.integer(substr(cellinfo[,1],2,1000))),];
	cellmap = match(cells, cellinfo[,1]);
	nn = dim(Y)[2];

	library("data.table");
	#egene = as.matrix(fread("/gpfs/group/yzz2/default/scratch/roadmap_analysis/gtex/gene_id.txt",header=F));
	#egene = as.matrix(fread("/gpfs/group/yzz2/default/scratch/GWAS/GSK/ToolVariant_cred-sets/gene_id.txt",header=F));
	egene = as.matrix(fread("/gpfs/group/yzz2/default/scratch/roadmap_analysis/PCHiC/gene_id.txt",header=F));
	l=as.integer(regexpr("\\.",egene[,1]));
	if(length(which(l>0))>0)
	{	egene[which(l>0)] = substr(egene[which(l>0),1],1,l[which(l>0)]-1);
	} else
	{	egene = egene[,1];
	}
	#snp = as.matrix(fread("/gpfs/group/yzz2/default/scratch/roadmap_analysis/gtex/var_id.txt",header=F));
	#snp = as.matrix(fread("/gpfs/group/yzz2/default/scratch/GWAS/GSK/ToolVariant_cred-sets/var_id.txt",header=F));
	snp = as.matrix(fread("/gpfs/group/yzz2/default/scratch/roadmap_analysis/PCHiC/var_id.txt",header=F));
	snpchr = as.matrix(snp[,1]);
	snppos = round(as.integer(as.matrix(snp[,2]))/200);	
	#pair = as.matrix(fread("/gpfs/group/yzz2/default/scratch/roadmap_analysis/gtex/allsignif_snpgene_pairs.txt"));
	#pair = as.matrix(fread("/gpfs/group/yzz2/default/scratch/GWAS/GSK/ToolVariant_cred-sets/allsignif_snpgene_pairs.txt"));
	pair = as.matrix(fread("/gpfs/group/yzz2/default/scratch/roadmap_analysis/PCHiC/allsignif_snpgene_pairs.txt"));
	#pair[,4] = -log10(pair[,4]);

	tt = which(chr == paste("chr",snpchr[pair[,2]],sep=""));
	pair = pair[tt,];

	for(i in chr)
	{	print(i);
		g = fread(paste(path,statepref,i,".state",sep=""));
		a=system(paste("wc -l ",path, statepref,i,".para",sep=""),intern=T);
		gn=as.integer(substr(a,1,regexpr(" ",a)[1]))-1;
		pos = as.matrix(g[,3]);
		state = as.matrix(g[,5:131])[,cellmap];
		t = which(gene[,2]==substr(i,4,1000));
		emap = match(snppos, pos / 200);

		print("state file loaded");
	
		gid = tinv = rinv = NULL;	
		for(j in t)
		{	tt = which(pos >= tss[j] - 1000 & pos < tss[j] + 1000);
			if(length(tt) == 0) next;
			rr = which(pos >= tss[j] - maxdist & pos < tss[j] + maxdist);
			if(length(rr) < 2) next;
			gid = c(gid, j);
			tinv = rbind(tinv, range(tt));
			rinv = rbind(rinv, range(rr));	
		}

		if(length(gid) < 10) next;
		#e = c(as.matrix(lm(c(Y[gid,])~as.factor(state[as.integer((tinv[,1]+tinv[,2])/2),])-1)$coef));
		#e[is.na(e)==T]=0;
		ee=NULL;
		str = paste(path,statepref,"coef",sep="");
		if(file.exists(str) == F)
		{	ee = getStateCoef(paste(path,statepref,sep=""));
		} else
		{	ee = read.table(str);
		}
		
		lab = rownames(ee);
		l = apply(t(t(lab)),1,function(x){tail(unlist(gregexpr("_",x)),n=1)});
		lab = substr(lab,1,l-1);
		ee = t(array(ee[,(dim(ee)[2]+1)/2],dim=c(gn,dim(ee)[1]/gn)));
		rownames(ee) = unique(lab);
		colnames(ee) = 1:gn-1;
		t=which(apply(ee,1,sd) == 0);
		if(length(t)>0) ee[t,] = ee[dim(ee)[1],];

		fout = paste("/gpfs/scratch/yzz2/rna/", fpref, ".", statepref,i,".r",sep="");
		gst = 0;
		for(j in (gst+1):length(gid))
		{	if(length(ee)>0) e = ee[which(gene[gid[j],6]==rownames(ee)),];			
			X = t(array(e[c(state[rinv[j,1]:rinv[j,2],])+1], dim=c(rinv[j,2]-rinv[j,1]+1,nn)));
			#X = smoothMatrix(X);
			yy = Y[gid[j],];
			r = cor(yy, X+rnorm(length(X),0,0.0001));
			r[is.na(r)==T] = 0;

			eq = rep(0,rinv[j,2]-rinv[j,1]+1+50);
			tt=which(gene[gid[j],1]==egene);
			if(length(tt)>0)
			{	tt = which(pair[,3]==tt[1]);
				tt = tt[which(is.na(emap[pair[tt,2]])==F)];
				if(length(tt)>0)
				{	tt=tt[which(emap[pair[tt,2]] >= rinv[j,1]-25 & emap[pair[tt,2]] <= rinv[j,2]+25)];
					if(length(tt)>0)
					{	eq[emap[pair[tt,2]]-rinv[j,1]+1+25] = pair[tt,4];
					}
				}
			}
			teq=NULL;for(i in -25:25) teq=cbind(teq,eq[1:(rinv[j,2]-rinv[j,1]+1)+25+i]);
			eq = teq[,26];
			eq5 = apply(teq[,26+-5:5],1,max);
			eq25 = apply(teq[,26+-25:25],1,max);
			
			print(j);
			#a = 1:dim(X)[2];
			a = which(r > minr);
			if(length(a) == 0) next;
			write.table(array(c(rep(as.matrix(gene[gid[j],1]),length(a)),(rinv[j,1]:rinv[j,2])[a],round(r[a]*1000)/1000, round(apply(X,2,mean)[a]*100)/100,round(apply(X,2,sd)[a]*100)/100, round(eq[a]*100)/100, round(eq5[a]*100)/100, round(eq25[a]*100)/100),dim=c(length(a),8)), fout, quote=F,row.names=F,col.names=F,append=(j>gst+1));
		}
	}
}

checkeqtl<-function(fcase, rmtss=FALSE, usemcmc=FALSE, subset=NULL, colid=6, suff="eqtl")
{	x=read.table(paste(fcase,".r",sep=""));
	#t=which(x[,colid]>0);x[t,colid]=10^(-x[t,colid]);

	if(!usemcmc)
	{	u=read.table(fcase,sep=" ");
	} else
	{	u=read.table(paste(fcase,".mcmc",sep=""));
	}
	if(length(subset)>0) u=u[subset[which(subset<=dim(u)[1])],];
	gene = read.table("Ensembl_v65.Gencode_v10.ENSG.gene_info");
	gtype = unique(as.matrix(gene[,6]));
	u = u[which(is.na(match(u[,1],as.matrix(gene[,1])))==F),];
	gt = gene[match(u[,1],as.matrix(gene[,1])),6];
	gtid = match(gt, gtype);

	tm=match(as.matrix(x[,1]),as.matrix(u[,1]));
	gtid = gtid[tm];
	mm=rep(0,dim(x)[1]);
	mask=rep(0,dim(x)[1]);
	for(i in 1:dim(u)[1])
	{	t=which(tm==i);
		if(length(t)==0) next;
		if(max(x[t,colid])<1e-10) mask[t]=1;
		if(!usemcmc)
		{	ll=as.integer(unlist(strsplit(as.matrix(u[i,9]),",")));
			pp=rep(1,length(ll));#as.numeric(unlist(strsplit(as.matrix(u[i,10]),",")));
			
		} else
		{	ll=as.integer(unlist(strsplit(as.matrix(u[i,2]),",")));
			pp=as.numeric(unlist(strsplit(as.matrix(u[i,3]),",")));
		}
		tt = t[match(ll,x[t,2])];
		ttt=which(is.na(tt)==F);
		if(length(ttt)>0)
		{	mm[tt[ttt]] = pp[ttt]; }
		if(i%%100==0) print(i);
	}
	if(rmtss)
	{	alltss=as.integer(read.table(paste("tsspos.chr",u[1,2],sep=""))[,3]);
		mask = (mask | alltss[x[,2]]);
	}

	a=a0=b=b0=n=n0=NULL;
	mr = as.integer((min(x[,3])+1)*10);
	mr = 11;
	for(i in mr:19)
	{	k=(i-10)/10;
		tt=which(x[,3]>k & x[,3]<=k+0.1 & mask==0);
		n0=rbind(n0,as.integer(tabulate(gtid[tt],nbins=length(gtype))));
		ta = tb = NULL;
		for(j in 1:length(gtype))
		{	ttt = which(gtid[tt]==j);
			ta[j] = sum(x[tt[ttt],colid]);
			tb[j] = sum(as.integer(x[tt[ttt],colid]>0));
		}
		ta[is.na(ta)==T]=0;
		tb[is.na(tb)==T]=0;
		a0=rbind(a0, ta);
		b0=rbind(b0, tb);

		tt=which(x[,3]>k & x[,3]<=k+0.1 & mask==0 & mm != 0);
		ta = tb = tn = NULL;
		for(j in 1:length(gtype))
		{	ttt = which(gtid[tt]==j);
			tn[j] = sum(mm[tt[ttt]]);
			ta[j] = sum(x[tt[ttt],colid]*mm[tt[ttt]]);
			tb[j] = sum(as.integer(x[tt[ttt],colid]>0)*mm[tt[ttt]]);
		}
		tn[is.na(tn)==T]=0
		ta[is.na(ta)==T]=0;
		tb[is.na(tb)==T]=0;
		n=rbind(n, tn);
		a=rbind(a, ta);
		b=rbind(b, tb);
	}
	rownames(a)=rownames(b)=rownames(n)=((mr:19)-10)/10;
	rownames(a0)=rownames(b0)=rownames(n0)=((mr:19)-10)/10;
	colnames(a)=colnames(b)=colnames(n)=gtype;
	colnames(a0)=colnames(b0)=colnames(n0)=gtype;
	write.table(cbind(a,b,n),paste(fcase,".",suff,"check",sep=""),quote=F,row.names=F);
	write.table(cbind(a0,b0,n0),paste(fcase,".",suff,"check0",sep=""),quote=F,row.names=F);
}

checkeqtl2<-function(fpref, chr, maxdist = 200000, rmtss=FALSE, usemcmc=FALSE, subset=NULL, colid=6)
{	x=read.table(paste(fpref, ".", chr,".r",sep=""));

	if(!usemcmc)
	{	u=read.table(paste("/gpfs/group/yzz2/default/scratch/roadmap_analysis/rna/step3.3", ".", chr,sep=""),sep=" ");
	} else
	{	u=read.table(paste(fpref, ".", chr, ".mcmc",sep=""));
	}
	if(length(subset)>0) u=u[subset[which(subset<=dim(u)[1])],];
	gene = read.table("Ensembl_v65.Gencode_v10.ENSG.gene_info");
	gtype = unique(as.matrix(gene[,6]));
	u = u[which(is.na(match(u[,1],as.matrix(gene[,1])))==F),];
	gt = gene[match(u[,1],as.matrix(gene[,1])),6];
	gtid = match(gt, gtype);

	tm=match(as.matrix(x[,1]),as.matrix(u[,1]));
	gtid = gtid[tm];
	mm=rep(0,dim(x)[1]);
	mask=rep(0,dim(x)[1]);
	for(i in 1:dim(u)[1])
	{	t=which(tm==i);
		if(length(t)==0) next;
		if(max(x[t,colid])<1e-10) mask[t]=1;
		if(!usemcmc)
		{	ll=as.integer(unlist(strsplit(as.matrix(u[i,9]),",")));
			pp=rep(1,length(ll));#as.numeric(unlist(strsplit(as.matrix(u[i,10]),",")));
			
		} else
		{	ll=as.integer(unlist(strsplit(as.matrix(u[i,2]),",")));
			pp=as.numeric(unlist(strsplit(as.matrix(u[i,3]),",")));
		}
		tt = t[match(ll,x[t,2])];
		ttt=which(is.na(tt)==F);
		if(length(ttt)>0)
		{	mm[tt[ttt]] = pp[ttt]; }
		if(i%%100==0) print(i);
	}
	if(rmtss)
	{	alltss=as.integer(read.table(paste("tsspos.chr",u[1,2],sep=""))[,3]);
		mask = (mask | alltss[x[,2]]);
	}
	tss = gene[,3];
	tss[which(gene[,5]<0)] = gene[which(gene[,5]<0),4];
	tss = tss[match(x[,1],as.matrix(gene[,1]))];
	pos = read.table(paste("tsspos.",chr,sep=""));
	pos = pos[x[,2],2];
	dt = (pos - tss) * as.integer(as.matrix(gene[match(as.matrix(x[,1]),as.matrix(gene[,1])),5]));
	dt = round(dt/5000 + 0.5);#round(log10(abs(dt)+1) * sign(dt) * 10);
	hsz = maxdist / 5000;
	dbase=-hsz:hsz;
	#dt[which(dt>=0 & dt< 30)] = 30;
	#dt[which(dt> -30 & dt <0)] = -30;
	#dbase = c((-53):(-30),30:53);
	dm = match(dt,dbase);

	tcode=(gtid-1)*length(dbase)+dm;

	L = max(x[,2]);
	b=b0=n=n0=NULL;
	for(i in dbase) 
	{	l00=l0=l=rep(0,L);
		l00[x[which(dt==i),2]]=1;
		apply(as.matrix(x[which((dt)==i),c(2,colid)]),1,function(x){l[x[1]]<<-max(l[x[1]],x[2]);return(NULL)});
		apply(cbind(as.matrix(x[,2]),mm)[which((dt)==i),],1,function(x){l0[x[1]]<<-max(l0[x[1]],x[2]);return(NULL)});
		b=c(b,length(which(l>0 & l0>0 & l00==1)));
		n=c(n,(length(which(l0>0 & l00==1))));
		b0=c(b0,length(which(l>0 & l00==1)));
		n0=c(n0,length(which(l00==1)));
		print(i);   
	}
	return(cbind(b,n,b0,n0));
}

prepareWig<-function(fpref)
{	RT=RT1=NULL;
	gene=as.matrix(read.table("/gpfs/group/yzz2/default/RoadmapEpigenomics/RNA-Seq/Ensembl_v65.Gencode_v10.ENSG.gene_info"));
	options(scipen=1000);
	for(i in c(1:22,"X","Y"))
	{	pos=read.table(paste("/gpfs/group/yzz2/default/RoadmapEpigenomics/RNA-Seq/tsspos.chr",i,sep=""));
		score=rep(0,dim(pos)[1]);
		y=read.table(paste(fpref,".chr",i,".mcmc",sep=""));
		apply(as.matrix(y),1,function(x)
		{	l=as.integer(unlist(strsplit(x[2],",")));
			p=as.numeric(unlist(strsplit(x[3],",")));
			p=p-score[l];
			t=which(p>0);
			if(length(t)>0)
			{	score[l[t]]<<- score[l[t]]+p[t];
			}
			return(NULL);
		});
		score[which(score>1000)]=1000;
		t=which(score>0);
		rt=cbind(pos[t,1:2],pos[t,2]+200,score[t]);
		RT=rbind(RT,rt);

		loc=NULL;
		name=NULL;
		y1=as.matrix(read.table(paste(fpref,".chr",i,sep="")));
		tt=match(y1[,1],gene[,1]);
		ttt=which(is.na(tt)==F & is.na(gene[tt,7])==F);
		y1[ttt,1]=gene[tt[ttt],7];	
		apply(as.matrix(y1),1,function(x)
		{	l=as.integer(unlist(strsplit(x[9],",")));
			if(length(l)>0 & is.na(x[1])==F)
			{	loc <<- rbind(loc,pos[l,]);
				name <<- c(name, rep(x[1],length(l)));
			}
			return(NULL);
		});
		RT1=rbind(RT1,cbind(loc[,1:2],loc[,2]+200,name));
		print(i); 
	}
	
	RT=paste(RT[,1],RT[,2],RT[,3],RT[,4]);
	RT=c(paste("track type=wiggle_0 name=\"",fpref,"\" description=\"",fpref,"\"",sep=""),RT);
	write.table(RT,paste(fpref,".mcmc.wig",sep=""),quote=F,row.names=F,col.names=F);
	system(paste("/storage/home/yzz2/work/tools/wigToBigWig ", fpref,".mcmc.wig /storage/home/yzz2/work/tools/bedtools2/genomes/human.hg19.genome ",fpref,".mcmc.bw",sep=""));

	write.table(RT1,paste(fpref,".sel.bed",sep=""),quote=F,row.names=F,col.names=F);
	system(paste("/storage/home/yzz2/work/tools/sort-bed ",fpref,".sel.bed > ", fpref,".sel.bed1",sep=""));
	system(paste("/storage/home/yzz2/work/tools/bedToBigBed ", fpref, ".sel.bed1 /storage/home/yzz2/work/tools/bedtools2/genomes/human.hg19.genome ",fpref,".sel.bb",sep=""));	
}

smoothMatrix<-function(X,hsz=2)
{	nX = NULL;
	l = dim(X)[2];
	if(l==1) return(X);
	w = c(1/(hsz:1+1),1,1/(1:hsz+1));
	for(i in 1:l)
	{	t = -hsz:hsz;
		t = t[which(i+t>0 & i+t<=l)];
		if(length(t)==1)
		{	nX = cbind(nX, X[,i]);
		} else
		{	nX = cbind(nX, apply(t(t(X[,i+t])*w[hsz+1+t]),1,sum)/sum(w[hsz+1+t]));
		} 
	}
	return(nX);
}

stepMulti_inc<-function(x, y, penalty, sel=NULL, sz = 0)
{	n = length(y);
	x = array(x, dim=c(n, length(x)/n));
	l = dim(x)[2];	
	#y = y - mean(y);
	#x = t(t(x) - apply(x, 2, mean));

	lp = getlp_multi(x,y,sel,n);
	
	for(i in 1:l)
	{	nlp = rep(lp - 1000000, l);
		t = 1:l;
		if(length(sel)>0) t=t[-sel];
		for(j in t)
		{	nlp[j] = getlp_multi(x,y,c(sel,j),n);
		}
		if(sz > 0 | max(nlp) > lp + penalty)
		{	k = which(nlp == max(nlp))[1];
			sel = c(sel, k);
			lp = nlp[k];
			cat(k,"");
		} else
		{	break;
		}
		if(sz > 0 & i >= sz) break;
	}
	return(sel);
}

stepMulti_dec<-function(x, y, penalty, sel=NULL)
{	n = length(y);
	x = array(x, dim=c(n, length(x)/n));
	l = dim(x)[2];	
	if(length(sel)==0) sel = 1:l;
	#y = y - mean(y);
	#x = t(t(x) - apply(x, 2, mean));

	lp = getlp_multi(x,y,sel,n);
	
	while(length(sel) > 0)
	{	nlp = rep(lp - 1000000, length(sel));
		for(j in 1:length(sel))
		{	nlp[j] = getlp_multi(x,y,sel[-j],n);
		}
		if(lp - penalty < max(nlp))
		{	k = which(nlp == max(nlp))[1];
			cat(sel[k], "");#lp, nlp[k],"\n");
			sel = sel[-k];
			lp = nlp[k];
		} else
		{	break;
		}
	}
	return(sel);
}

checkOverlap<-function(fname, fnull=NULL, eqtl=FALSE, cutoff=5)
{	library("data.table");
	gene = as.matrix(read.table(fname, sep=" "));
	l = as.integer(regexpr("\\.",as.matrix(gene[,1])));
	if(length(which(l>0))>0)
	{	gene[which(l>0),1]=substr(as.matrix(gene[which(l>0),1]),1,l[which(l>0)]-1);
	}
	
	if(length(fnull)==0)
	{	fpref = paste(fname, ".bed", sep="");
		x = as.matrix(fread(fpref));
	} else
	{	fpref = fnull;
		x = as.matrix(fread(fnull));
	}
	
	if(eqtl)
	{	egene = as.matrix(fread("/gpfs/group/yzz2/default/scratch/roadmap_analysis/gtex/gene_id.txt",header=F));
	} else
	{	egene = as.matrix(fread("/gpfs/group/yzz2/default/scratch/roadmap_analysis/PCHiC/gene_id.txt",header=F));
	}
	l=as.integer(regexpr("\\.",egene[,1]));
	if(length(which(l>0))>0)
	{	egene[which(l>0)] = substr(egene[which(l>0),1],1,l[which(l>0)]-1);
	} else
	{	egene = egene[,1];
	}
	if(eqtl)
	{	snp = as.matrix(fread("/gpfs/group/yzz2/default/scratch/roadmap_analysis/gtex/var_id.txt",header=F));
		snppos = cbind(as.numeric(snp[,2]),as.numeric(snp[,2]));
	} else
	{	snp = as.matrix(fread("/gpfs/group/yzz2/default/scratch/roadmap_analysis/PCHiC/var_id.txt",header=F));
		snppos = cbind(as.numeric(snp[,2]),as.numeric(snp[,3]));
	}
	snpchr = as.matrix(snp[,1]);
	
	if(eqtl)
	{	pair = as.matrix(fread("/gpfs/group/yzz2/default/scratch/roadmap_analysis/gtex/allsignif_snpgene_pairs.txt"));
		pair[,4]=-log10(pair[,4]);
	} else
	{	pair = as.matrix(fread("/gpfs/group/yzz2/default/scratch/roadmap_analysis/PCHiC/allsignif_snpgene_pairs.txt"));
	}

	chr = as.matrix(gene[1,2]);
	gpos = as.matrix(read.table(paste("/storage/home/yzz2/work/RoadmapEpigenomics/data/chr",chr,".txt",sep=""))[,-1]);
	
	xpos = gpos[x[,3],1];

	dd = ss = rep(-1, dim(x)[1]);
	for(i in unique(x[,1]))
	{	t = which(x[,1]==i);
		te = which(egene == gene[i,1]);
		if(length(te)>0)
		{	ll = which(pair[,3]==te & pair[,4]>=cutoff);
			lll = length(ll);
			if(lll>0)
			{	v = array(snppos[pair[ll,2],],dim=c(lll,2));
				for(j in t)
				{	d = apply(v - xpos[j],1,function(z){if(sign(z[1])*sign(z[2])<=0){return(0)}else{return(min(abs(z)));}});
					id = which(d==min(d));
					ss[j] = max(pair[ll[id],4]);
					dd[j] = min(d);
				}
			}
		}
		print(i);
	}

	if(eqtl)
	{	write.table(cbind(dd,ss), paste(fpref,".eqtloverlap",sep=""),quote=F,row.names=F,col.names=F);
	} else
	{	write.table(cbind(dd,ss), paste(fpref,".pchicoverlap",sep=""),quote=F,row.names=F,col.names=F);
	}
}

summarizePCHiC<-function(fpref,fnull)
{
	X = X0 = Y = Y0 = NULL;
	for(i in paste("chr",c(1:22,"X"),sep=""))
	{	x=as.matrix(read.table(paste(fpref,".",i,".bed.pchicoverlap",sep="")));
		x0=as.matrix(read.table(paste(fnull,".",i,".bed.pchicoverlap",sep="")));
		y=as.matrix(read.table(paste(fpref,".",i,".bed",sep="")));
		y0=as.matrix(read.table(paste(fnull,".",i,".bed",sep="")));
		t=which(x[,1]>=0);
		t0=which(x0[,1]>=0);
		X=rbind(X,x[t,]);
		X0=rbind(X0,x0[t0,]);
		Y=rbind(Y,y[t,]);
		Y0=rbind(Y0,y0[t0,]);
		print(i);
	}

	A = A0 = N = N0 = NULL;
	for(j in 1:25)
	{	a=a0=n=n0=NULL;
		for(i in 1:100) 
		{	a[i]=(length(which(X[,1]<=j*200 & abs(Y[,4])>i*10000)));
			n[i]=(length(which(abs(Y[,4])>i*10000)));
			a0[i]=(length(which(X0[,1]<=j*200 & abs(Y0[,4])>i*10000)));
			n0[i]=(length(which(abs(Y0[,4])>i*10000)));
		}
		A = rbind(A, a);
		N = rbind(N, n);
		A0 = rbind(A0, a0);
		N0 = rbind(N0, n0);	
		print(j);
	}
	
	rt = NULL;
	rt$A = A;
	rt$N = N;
	rt$A0 = A0;
	rt$N0 = N0;
	return(rt);
}

getStateEffect<-function(fpref, chr, shift=0)
{
	options(scipen=1000);
	gene=read.csv("Ensembl_v65.Gencode_v10.ENSG.gene_info",sep="\t",header=F);
	Y=array(0,dim=c(dim(gene)[1],56));
	x=read.table("57epigenomes.RPKM.pc",header=T);
	cells=colnames(x)[-(1:2)];
	m=match(gene[,1],x[,1]);
	t=which(is.na(m)==F);
	Y[t,]=as.matrix(x[m[t],-(1:2)]);
	x=read.table("57epigenomes.RPKM.nc",header=T);
	m=match(gene[,1],x[,1]);
	t=which(is.na(m)==F);
	Y[t,]=as.matrix(x[m[t],-(1:2)]);
	x=read.table("57epigenomes.RPKM.rb",header=T);
	m=match(gene[,1],x[,1]);
	t=which(is.na(m)==F);
	Y[t,]=as.matrix(x[m[t],-(1:2)]);

	Y=log2(Y+0.01);

	t=which(apply(Y,1,sd) > 0);
	Y=Y[t,];
	gene=gene[t,];	
	genechr = paste("chr",gene[,2],sep="");
	tss = gene[,3];
	gsign = gene[,5];
	tss[gene[,5]<0] = gene[gene[,5]<0,4];

	cellinfo = read.table("/gpfs/group/yzz2/default/RoadmapEpigenomics/extra/epigenomeID.txt",header=T,comment="!",sep="\t");
	cellinfo = cellinfo[order(as.integer(substr(cellinfo[,1],2,1000))),];
	cellmap = match(cells, cellinfo[,1]);
	nn = dim(Y)[2];

	library("data.table");
	C = P = G = NULL;
	for(i in chr)
	{	print(i);
		t = which(genechr == i);
		if(length(t)==0) next;

		g = fread(paste(fpref,i,".state",sep=""));
		tchr = as.matrix(g[,2]);
		pos = as.matrix(g[,3]);
		g = as.matrix(g[,5:131]);
		g = g[,cellmap];
		G = rbind(G, g);
		C = c(C, tchr);
		P = c(P, pos);
	}	
	gn = max(G)+1;

	rt = NULL;
	for(j in shift)
	{	tG = tP = NULL;
		for(i in chr)
		{	print(c(j,i));
			t = which(C == i);
			tg = which(genechr == i);
			if(length(tg) == 0 | length(t) == 0) next;
			tp = 200 * (as.integer(tss[tg] / 200) + j * gsign[tg]);
			m = match(tp, P[t]);
			tt = which(is.na(m)==F);
			if(length(tt) == 0) next;
			tG = c(tG, tg[tt]);
			tP = c(tP, t[m[tt]]);
		}
		rt$Y = rbind(rt$Y, Y[tG,]);
		rt$X = rbind(rt$X, G[tP,]);
		rt$G = rbind(rt$G, gene[tG,]);
		rt$shift = c(rt$shift, rep(j,length(tG)));
	}
	return(rt);
}

runStateEffect<-function(fpref, shift=0)
{	rt = getStateEffect(fpref, paste("chr",c(1:22,"X","Y"),sep=""), shift);
	gn = max(rt$X)+1;
	gname = unique(as.matrix(rt$G[,6]));

	rr = ee = NULL;
	for(i in shift)
	{	trr = tee = NULL;
		for(j in gname)
		{	print(c(i,j));
			t = which(rt$shift == i & rt$G[,6]==j);
			#if(length(t)>0) rt$Y[t,] = rt$Y[t,] - mean(c(rt$Y[t,]));
			if(length(t)< 100)
			{	trr = c(trr, 0);
				tee = c(tee, rep(0,gn));
				next;
			}
			tr=lm(c(rt$Y[t,])~as.factor(c(rt$X[t,]))-1);	
			trr = c(trr,summary(tr)$adj.r.square);
			a = sort(unique(c(rt$X[t,])));
			ttt = rep(0,gn);ttt[a+1] = coef(tr);
			tee = c(tee,ttt);
		}
		t = which(rt$shift == i);
		tr=lm(c(rt$Y[t,])~as.factor(c(rt$X[t,]))-1);	
		trr = c(trr,summary(tr)$adj.r.square);
		a = sort(unique(c(rt$X[t,])));
		ttt = rep(0,gn);ttt[a+1] = coef(tr);
		tee = c(tee,ttt);

		rr = cbind(rr, trr);
		ee = cbind(ee, tee);
	}
	colnames(rr) = colnames(ee) = shift;
	rownames(rr) = c(gname, "all");
	rownames(ee) = paste(rep(c(gname, "all"), each=gn), rep(1:gn-1,length(gname)+1),sep="_");
	
	rt$rr = rr;
	rt$ee = ee;
	return(rt);
}

getStateCoef<-function(statepref)
{	gene=read.table("Ensembl_v65.Gencode_v10.ENSG.gene_info");
	Y=array(0,dim=c(dim(gene)[1],56));
	x=read.table("57epigenomes.RPKM.pc",header=T);
	cells=colnames(x)[-(1:2)];
	m=match(gene[,1],x[,1]);
	t=which(is.na(m)==F);
	Y[t,]=as.matrix(x[m[t],-(1:2)]);
	x=read.table("57epigenomes.RPKM.nc",header=T);
	m=match(gene[,1],x[,1]);
	t=which(is.na(m)==F);
	Y[t,]=as.matrix(x[m[t],-(1:2)]);
	x=read.table("57epigenomes.RPKM.rb",header=T);
	m=match(gene[,1],x[,1]);
	t=which(is.na(m)==F);
	Y[t,]=as.matrix(x[m[t],-(1:2)]);

	Y=log2(Y+0.01);

	t=which(apply(Y,1,sd) > 0);
	Y=Y[t,];
	gene=gene[t,];	
	genechr = paste("chr",gene[,2],sep="");
	tss = gene[,3];
	tss[gene[,5]<0] = gene[gene[,5]<0,4];

	library("MASS");
	cellinfo = read.table("/gpfs/group/yzz2/default/RoadmapEpigenomics/extra/epigenomeID.txt",header=T,comment="!",sep="\t");
	cellinfo = cellinfo[order(as.integer(substr(cellinfo[,1],2,1000))),];
	cellmap = match(cells, cellinfo[,1]);
	celln = dim(Y)[2];
	library("data.table");
	a=system(paste("wc -l ",statepref,"chr1.para",sep=""),intern=T);
	gn=as.integer(substr(a,1,regexpr(" ",a)[1]))-1;
	X = array(0, dim=c(dim(Y)[1]*celln,gn));
	sel = rep(1,dim(Y)[1]);
	for(i in paste("chr",c(1:22,"X","Y"),sep=""))
	{	print(i);
		Z = NULL;
		g=fread(paste(statepref,i,".state",sep=""));
		pos = as.matrix(g[,3]);
		state = as.matrix(g[,5:131])[,cellmap];

		t = which(genechr == i);
		for(j in t)
		{	tt = which(pos >= tss[j] - 1000 & pos < tss[j] + 1000);
			if(length(tt) > 0)
			{	X[(j-1)*celln+1:celln,] = t(apply(array(state[tt,],dim=c(length(tt),celln)),2,function(z){tabulate(z+1,nbins=gn)})) / length(tt);
				
			} else { sel[j] = 0;	}
			print(c(which(j==t),length(t)));
		}
	}

	geneid = unique(gene[,6]);
	ee = array(0, dim=c(length(geneid),gn));
	for(i in 1:length(geneid))
	{	t = which(gene[,6] == geneid[i] & sel == 1);
		if(length(t)>0)
		{	ty = c(t(Y[t,]));
			tx = X[rep((t-1)*celln,each=celln) + rep(1:celln,length(t)),];
			e = c(as.matrix(lm(ty~tx-1)$coef));
			e[is.na(e)==T]=0;
			ee[i,] = e;
		}
	}
	ee = array(t(ee),dim=c(length(ee),1));
	rownames(ee) = paste(rep(as.matrix(geneid),each=gn),rep(1:gn-1,length(geneid)),sep="_");
	write.table(round(ee*10000)/10000, paste(statepref,"coef",sep=""),quote=F);
	return(ee);
}
