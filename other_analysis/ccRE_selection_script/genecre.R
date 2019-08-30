genereg<-function(rna, state, pair, tssdist, prior=NULL)
{
	celln = dim(rna)[2];
	staten = dim(state)[2]/celln;
#	if(length(prior) == 0) 
#	{	prior = array(1, dim=rep(max(c(pair[,1:2])),2));
#	}
	l = dim(prior)[1];
	y = x = x0 = NULL;
	ut = unique(pair[,3]);
	for(i in 1:length(ut))
	{	tt = which(pair[,3] == ut[i]);
		if(length(tt)==0) next;
		y = rbind(y, rna[ut[i],]);
		ttt = pair[tt[1],1]+-tssdist:tssdist;
		ttt = ttt[ttt>0 & ttt <= dim(state)[1]];
		x0 = rbind(x0, apply(array(state[ttt,],dim=c(length(ttt),dim(state)[2])),2,mean));
		if(length(prior) > 0)
		{	x = rbind(x, apply(array(state[pair[tt,2],]*prior[(pair[tt,1]-1)*l+pair[tt,2]],dim=c(length(tt),staten*celln)),2,mean));
		} else
		{	x = rbind(x, apply(array(state[pair[tt,2],],dim=c(length(tt),staten*celln)),2,mean));
		}
		cat(i,",",sep="");
	}
	xx = xx0 = NULL;
	for(i in 1:celln)
	{	xx = rbind(xx, array(x[,(i-1)*staten + 1:staten],dim=c(dim(x)[1],staten)));
		xx0 = rbind(xx0, array(x0[,(i-1)*staten + 1:staten],dim=c(dim(x0)[1],staten))); 
	}
	print("");

	rt = NULL;
	rt$y = c(y);
	rt$z = c((y-apply(y,1,mean))/(apply(y,1,sd)+1e-3));
	rt$x = xx;
	rt$x0 = xx0;
	return(rt);
}

refineList<-function(y,x,x0,state,pair,sel,lessone, itern)
{	ut = unique(pair[,3]);
	n = length(ut);
	celln = length(y)/n;
	staten = dim(x)[2];
	tx = cbind(log2(x+0.001),log2(x0+0.001));
	e = f = NULL;
	for(i in (1:celln)[-lessone])
	{	tt = rep((c(i, lessone)-1) * n, each=n) + rep(1:n,2);
		r = lm(y[-tt]~tx[-tt,]);
		te = r$coef;
		te[is.na(te)==T] = 0;
		e = cbind(e, te);
		f = c(f, c(cbind(1,tx[(i-1)*n + 1:n,])%*%te));
	}	
	mm = c(apply(array((y[-((lessone-1)*n+1:n)]-f)^2,dim=c(n, celln - 1)),1,sum));
	nx = x;

	for(i in 1:n)
	{	t = which(pair[,3] == ut[i]);
		me=NULL;
		if(length(t) > 0)
		{	for(j in 1:length(t))
			{	ttt = (array(x[(1:celln-1)*n+i,],dim=c(celln,staten)) * sum(sel[t]) - (2 * sel[t[j]] - 1) * t(array(state[pair[t[j],2],],dim=c(staten, celln)))) / (sum(sel[t]) - (2 * sel[t[j]] - 1) + 1e-10);
				ttt[ttt<0] = 0;
				f = cbind(1,log2(ttt+0.001),log2(x0[(1:celln-1)*n+i,]+0.001))%*%e;
				f = diag(f[-lessone,]);
				me[j] = sum((y[((1:celln)[-lessone]-1)*n+i]-f)^2);
			}
			#j = which(me == min(me))[1];
			#if(me[j] < mm[i])
			j = which(me < mm[i]);
#j=rbinom(length(me),1,prob=exp(mm[i]-me)/(1+exp(mm[i]-me)));
#j=which(j==1);
			if(length(j)>0)
			{	#ttt = (x[(1:celln-1)*n+i,] * sum(sel[t]) - (2 * sel[t[j]] - 1) * t(array(state[pair[t[j],2],],dim=c(staten, celln)))) / (sum(sel[t]) - (2 * sel[t[j]] - 1) + 1e-10);
				#mm[i] = me[j];

				tp = exp((mm[i]-me[j]-max(mm[i]-me[j]))/2);
				tp[is.na(tp)]=0;
				tp = tp + 1e-10/length(j);
				if(length(j)>round(1000/itern)+1) j = j[sample(length(j),size=round(100/itern)+1,prob=tp)];
				ttt = (array(x[(1:celln-1)*n+i,],dim=c(celln,staten)) * sum(sel[t]) - t(array(apply(array(state[pair[t[j],2],]*(2*sel[t[j]]-1),dim=c(length(j),dim(state)[2])),2,sum),dim=c(staten, celln)))) / (sum(sel[t]) - sum(2 * sel[t[j]] - 1)+1e-10);
				ttt[ttt<0] = 0;
				f = cbind(1, log2(ttt+0.001),log2(x0[(1:celln-1)*n+i,]+0.001))%*%e;
				f = diag(f[-lessone,]);
				mm[i] = sum((y[((1:celln)[-lessone]-1)*n+i]-f)^2);

				nx[(1:celln-1)*n+i,] = ttt;
				sel[t[j]] = 1 - sel[t[j]];
			}
		}
	}
	rt = NULL;
	rt$x = nx;
	rt$sel = sel;
	return(rt);
}

runRefine<-function(rt, celln, ss, pair, lessone, B = 100)
{	r=NULL;
	r$x=rt$x;
	k = length(rt$y) / celln;
	a00=summary(lm(rt$y[(lessone-1)*k+1:k]~log2(rt$x0[(lessone-1)*k+1:k,]+0.001)))$adj.r.squared;
	ma=summary(lm(rt$y[(lessone-1)*k+1:k]~log2(r$x[(lessone-1)*k+1:k,]+0.001)+log2(rt$x0[(lessone-1)*k+1:k,]+0.001)))$adj.r.squared;
	ms=sel=rep(1,dim(pair)[1]);
	a0 = a = ma;n0 = n = mean(sel);
	print(c(0,n,a));
	for(i in 1:B)
	{	r = refineList(rt$y,r$x,rt$x0,ss,pair, sel, lessone, i);
		if(length(which(sel != r$sel)) == 0) break;
		sel = r$sel;
		a[i]=summary(lm(rt$y[(lessone-1)*k+1:k]~log2(r$x[(lessone-1)*k+1:k,]+0.001)+log2(rt$x0[(lessone-1)*k+1:k,]+0.001)))$adj.r.squared;
		n[i]=mean(sel);
		print(c(i,n[i],a[i]));
		if(a[i] > ma)
		{	ma = a[i];
			ms = sel;
		}
	}
	r$n = n;
	r$n0 = n0;
	r$a = a;
	r$a0 = a0;
	r$a00 = a00;
	r$ma = ma;
	r$msel = ms;
	r$sel = sel;
	return(r);
}

run<-function(rna, tss, tssdist, statepref, chr, e = NULL, lessone = 12, maxdist=1e6, cut=0.2, B=100, fixEffect=FALSE)
{	print("load data ...");
	tss = as.integer(tss/200)+1;
	library("data.table")
	state = fread(paste(statepref,".state",sep=""));
	scr = as.matrix(state[,2]);
	pos = as.integer(as.matrix(state[,3])/200)+1;
	state = as.matrix(state[,-(1:4)]);
	state = state[,match(colnames(rna),colnames(state))];
	celln = dim(state)[2];
	t = which(scr == chr);
	pos = pos[t];
	state = state[t,];

	G = max(state)+1;
	t = tabulate(c(state)+1,nbins=G);
	k = which(t == max(t))[1] - 1;
	l = max(c(tss, pos));
	ss = array(k, dim=c(l,celln));
	ss[pos,] = state;

	tcre = read.table("vision_cres.txt");
	tcre = as.matrix(tcre[which(tcre[,1]==chr),2:3]);
	cre = rep(0, l);
	if(length(tcre) > 2)
	{	tcre[,1] = as.integer(tcre[,1]/200)+1;
		tcre[,2] = as.integer((tcre[,2]+199)/200)+1;
		for(i in 1:dim(tcre)[1])
		{	cre[max(1,tcre[i,1]):min(l,tcre[i,2])] = i;
			#cre[as.integer((tcre[i,1]+tcre[i,2])/2)] = 1;
		}
	}

	print("find gene-cre pairs ...");
	if(length(e) == 0)
	{	sss = array(0, dim=c(length(tss) * celln,G));
		for(i in 1:length(tss))
		{	ttt = tss[i] + -tssdist:tssdist;
			ttt = ttt[ttt > 0 & ttt <= dim(ss)[1]];
			aaa = t(apply(array(ss[ttt,],dim=c(length(ttt),dim(ss)[2])),2,function(z){tabulate(z+1,G)/length(z)}));
			sss[(1:celln-1)*length(tss)+rep(i,celln),] = aaa;
		}
		e = lm(c(rna)~sss-1)$coef;
		e[which(is.na(e)==T)] = 0;
	} 
	sse = array(e[ss+1],dim=dim(ss));

#	if(length(tcre) > 2)
#	{	osse = sse;
#		for(i in 1:dim(tcre)[1])
#		{	sse[as.integer((tcre[i,1]+tcre[i,2])/2),] = apply(array(osse[tcre[i,1]:tcre[i,2],],dim=c(tcre[i,2]-tcre[i,1]+1,dim(sse)[2])),2,mean);
#		}
#	}

	tr = (rna - apply(rna, 1, mean)) / (apply(rna,1,sd)+1e-5) / sqrt(celln-1);
	ts = (sse - apply(sse, 1, mean)) / (apply(sse,1,sd)+1e-5) / sqrt(celln-1);
	pair = NULL;
	for(i in 1:length(tss))
	{	a = max(1, tss[i]-maxdist/200):min(l, tss[i]+maxdist/200);
		rr = c(tr[i,]%*%t(ts[a,]));
		trr = c(tr[i,-lessone]%*%t(ts[a,-lessone]));
		t = which(trr >=cut & cre[a]>0);
		if(length(t)>0)
		{	nt = NULL;
			for(j in unique(cre[a[t]]))
			{	tt = t[which(cre[a[t]] == j)];
				nt = c(nt, tt[which(trr[tt]==max(trr[tt]))[1]]);
			}
			t = nt;
		}
		t = unique(c(which(a == tss[i]),t));#which(trr>=cut & cre[a]>0)));
		pair = rbind(pair, cbind(tss[i], a[t], i, rr[t]));
	}	
	print(dim(pair));

	print("prepare training ...");
	kk = rep(0, G);
	kk[k + 1] = 1;
	ss = array(rep(rep(kk, each=l),celln), dim=c(l, celln * G));
	for(i in 1:dim(state)[1])
	{	tt = rep(0, celln * G);
		tt[(1:celln-1)*G+state[i,]+1] = 1;
		ss[pos[i],] = tt;
	}
	#if(length(tcre) > 2)
	#{	oss = ss;
	#	for(i in 1:dim(tcre)[1])
	#	{	ss[as.integer((tcre[i,1]+tcre[i,2])/2),] = apply(array(oss[tcre[i,1]:tcre[i,2],],dim=c(tcre[i,2]-tcre[i,1]+1,dim(sse)[2])),2,mean);
	#	}
	#}
if(fixEffect) 
{ ss=2^sse; }
	rt = genereg(rna, ss, pair, tssdist, prior=NULL);

	print("select cres ...");
	r = runRefine(rt, celln, ss, pair, lessone, B = B);
	rt$nx = r$x;
	rt$sel = r$sel;
	rt$msel = r$msel;
	rt$n = r$n;
	rt$n0 = r$n0;
	rt$a = r$a;
	rt$a0 = r$a0;
	rt$a00 = r$a00;
	rt$ma = r$ma;
	rt$pair = pair;

	return(rt);	
}

sumCres<-function(fpref)
{	rrt = NULL;
	rrt$N0 = rrt$NN = rrt$G = rrt$A0 = rrt$AA = rrt$Ln = rrt$L1 = rrt$chr = NULL;
	for(k in 1:4)
	{	print(k);
		N0 = N = G = A0 = A = NULL;
		for(i in paste("chr",c(1:19,"X"),sep=""))
		{	n0 = n = g = a0 = a = NULL;
			for(j in 1:12)
			{	load(paste(fpref,".",i,".",k,".",j,".Rdat",sep=""));
				n0 = c(n0, dim(rt$pair)[1]);
				n = c(n, sum(rt$msel));
				g = c(g, length(unique(rt$pair[,1])));
				a0 = c(a0, rt$a0);
				a = c(a, rt$ma);
				if(length(which(rrt$chr == i))==0) 
				{	rrt$chr = c(rrt$chr, rep(i, max(rt$pair[,2])*1.1));
					rrt$Ln = rbind(rrt$Ln, array(0,dim=c(max(rt$pair[,2])*1.1,4)));
					rrt$L1 = rbind(rrt$L1, array(0,dim=c(max(rt$pair[,2])*1.1,4)));
				}
				t = which(rrt$chr == i);
				tt = table(rt$pair[rt$msel==1,2]);
				rrt$Ln[t[as.integer(names(tt))],k] = rrt$Ln[t[as.integer(names(tt))],k] + as.integer(tt);
				rrt$L1[t[as.integer(names(tt))],k] = rrt$L1[t[as.integer(names(tt))],k] + 1;	 
			}
			N0 = rbind(N0, n0);
			N = rbind(N, n);
			G = rbind(G, g);
			A0 = rbind(A0, a0);
			A = rbind(A, a);
		}
		rrt$N0 = cbind(rrt$N0, N0);
		rrt$NN = cbind(rrt$NN, N);
		rrt$G = cbind(rrt$G, G);
		rrt$A0 = cbind(rrt$A0, A0);
		rrt$AA = cbind(rrt$AA, A);
	}		
	nchr = npos = nL = nL1 = NULL;
	for(i in paste("chr",c(1:19,"X"),sep=""))
	{	print(i);
		t = which(rrt$chr == i);
		tt = which(apply(rrt$Ln[t,],1,sum) > 0);
		nchr = c(nchr, rrt$chr[t[tt]]);	
		npos = c(npos, (1:length(t))[tt]);
		nL = rbind(nL, rrt$Ln[t[tt],]);
		nL1 = rbind(nL1, rrt$L1[t[tt],]);
	}
	rrt$chr = nchr;
	rrt$pos = npos;
	rrt$Ln = nL;
	rrt$L1 = nL1;

	cres=read.table("vision_cres.txt");
	cpos=cbind(as.integer(cres[,2]/200)+1,as.integer(cres[,3]/200)+1);
	sel1=seln=array(0,dim=c(dim(cres)[1],4));
	for(i in paste("chr",c(1:19,"X"),sep=""))
	{	t=which(cres[,1]==i);
		tt=which(rrt$chr==i);
		l1=ln=array(0,dim=c(max(c(rrt$pos[tt],cpos[t,2])),4));
		l1[rrt$pos[tt],]=rrt$L1[tt,];
		ln[rrt$pos[tt],]=rrt$Ln[tt,];
		for(j in t)
		{	sel1[j,]=c(apply(array(l1[cpos[j,1]:cpos[j,2],],dim=c(cpos[j,2]-cpos[j,1]+1,4)),2,max));
			seln[j,]=c(apply(array(ln[cpos[j,1]:cpos[j,2],],dim=c(cpos[j,2]-cpos[j,1]+1,4)),2,max));
		}
		print(i);
	}
	rrt$sel1=sel1;
	rrt$seln=seln;
	return(rrt);
}

runatac<-function(rna, tss, tssdist, chr, lessone = 12, maxdist=1e6, cut=0.2, B=100)
{	print("load data ...");
	tss = as.integer(tss/200)+1;
	atac = read.table("vision_cres.mat.atacsig.txt",header=T,comment="!")
	achr = as.matrix(atac[,1]);
	pos = as.matrix(atac[,2:3]);
	t = which(achr == chr);
	pos = pos[t,];
	atac = as.matrix(atac[t,match(colnames(rna),colnames(atac))]);
	pos[,1] = as.integer(pos[,1]/200)+1;
	pos[,2] = as.integer((pos[,2]+199)/200)+1;
	celln = dim(atac)[2];
	l = max(c(tss, pos[,2]));
	ss = array(0, dim=c(l,celln));
	for(k in 0:max(pos[,2]-pos[,1]))
	{	t = which(pos[,1]+k <= pos[,2]);
		if(length(t) > 0)
			ss[pos[t,1]+k,] = 2^atac[t,];
	}

	tr = (rna - apply(rna, 1, mean)) / (apply(rna,1,sd)+1e-5) / sqrt(celln-1);
	ts = (ss - apply(ss, 1, mean)) / (apply(ss,1,sd)+1e-5) / sqrt(celln-1);
	pair = NULL;
	for(i in 1:length(tss))
	{	a = max(1, tss[i]-maxdist/200):min(l, tss[i]+maxdist/200);
		rr = c(tr[i,]%*%t(ts[a,]));
		trr = c(tr[i,-lessone]%*%t(ts[a,-lessone]));
		t = unique(c(which(a == tss[i]),which(trr>=cut)));
		pair = rbind(pair, cbind(tss[i], a[t], i, rr[t]));
	}	

	print("prepare training ...");
	rt = genereg(rna, ss, pair, tssdist, prior=NULL);

	print("select cres ...");
	r = runRefine(rt, celln, ss, pair, lessone, B = B);
	rt$nx = r$x;
	rt$sel = r$sel;
	rt$msel = r$msel;
	rt$n = r$n;
	rt$n0 = r$n0;
	rt$a = r$a;
	rt$a0 = r$a0;
	rt$a00 = r$a00;
	rt$ma = r$ma;
	rt$pair = pair;

	return(rt);	
}
