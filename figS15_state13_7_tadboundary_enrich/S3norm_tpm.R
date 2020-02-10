d = read.table('../rna/rnaTPM.txt', header=T)
dsig = d[,-c(1:4)]
dsig_linear = 2^dsig
dsig_linear[dsig_linear==dsig_linear[1,1]] = 0

median_mean = median(apply(dsig, 2, mean))
median_sd = median(apply(dsig, 2, sd))

get_snr = function(x, thresh){
x_fdr = get_fdr(x)
print(summary(x_fdr))
x_fdr_pk = x_fdr<thresh
x_fdr_bg = x_fdr>=thresh
snr = mean(x[x_fdr_pk]) / mean(x[x_fdr_bg])
return(snr)
}

get_frip = function(x, thresh){
x_fdr = get_fdr(x)
x_fdr_pk = x_fdr<thresh
frip = sum(x[x_fdr_pk]) / sum(x)
return(frip)
}


###### get NB model prob and size
get_true_NB_prob_size = function(x){
m=mean(x[x>0]);
m2=mean(x[x>0]^2);
p0 = length(which(x==0)) / length(x);
p = m/(m2-m^2 * (1-p0));
s = m * (1 - p0) * p /(1-p);
rt=c(p,s,p0);
for(i in 1:100){
op = p;
os = s;
p0=p^s;
#print(p0)
p=m/(m2-m^2*(1-p0));
if (p<0.001){
	p = 0.001
}
if (p>=0.999){
	p = 0.999
}
s=m * (1 - p0) * p / (1-p);
#rt=rbind(rt,c(p,s,p0));
rt = c(p,s,p0)
if(abs(op-p)<0.00001 & abs(os-s)<0.00001) break;
}
#print('change best_p0: ')
#print(p0)
return(rt);
}

get_pval = function(N, l, sig_0_size, sig_0_prob, num_0){
if (N != 0){
pval_new = pnbinom(N-1, sig_0_size, sig_0_prob, lower.tail=FALSE) / pnbinom(0, sig_0_size, sig_0_prob, lower.tail=FALSE) * (l-num_0)/l
} else {
pval_new = 1.0
}
return(pval_new)
}


get_fdr = function(sig){
### get NB
sig_notop = sig[sig<=quantile(sig, 0.99)]
obs_0_num = sum(sig_notop==0)
sig_notop_non0 = sig_notop[sig_notop>0]
sig_probT_sizeT = get_true_NB_prob_size(sig_notop_non0)
### get NB parameters
p0 = sig_probT_sizeT[3]
sig_size = sig_probT_sizeT[2]
sig_prob = sig_probT_sizeT[1]
### set limit for prob
if (sig_prob<0.001){
sig_prob = 0.001
}
if (sig_prob>=0.999){
sig_prob = 0.999
}
bin_num = length(sig)
### get NBP
sig = cbind(sig)
nb_pval = apply(sig, 1, function(x) get_pval(x, bin_num, sig_size, sig_prob, obs_0_num))
nb_pval[nb_pval<=1e-323] = 1e-323
nb_pval[nb_pval > 1] = 1
nb_pval[sig == 0] = 1
### get FDR NBP
xNBPfdr = p.adjust(nb_pval, 'fdr')
return(xNBPfdr)
}

s3norm = function(xref, xtar, thresh){
### ref
xref_fdr = get_fdr(xref)
xref_fdr_pk = xref_fdr<thresh
xref_fdr_bg = xref_fdr>=thresh
### tar
xtar_fdr = get_fdr(xtar)
xtar_fdr_pk = xtar_fdr<thresh
xtar_fdr_bg = xtar_fdr>=thresh
### cpk & cbg
cpk = (xref_fdr_pk * xtar_fdr_pk) !=0
cbg = (xref_fdr_bg * xtar_fdr_bg) !=0
print(sum(cpk))
print(sum(cbg))
### cpk cbg sig
xref_cpk = xref[cpk]
xtar_cpk = xtar[cpk]
xref_cbg = xref[cbg]
xtar_cbg = xtar[cbg]
### non0 mean
xref_cpk_non0_mean = mean(xref_cpk[xref_cpk>0])
xtar_cpk_non0_mean = mean(xtar_cpk[xtar_cpk>0])
xref_cbg_non0_mean = mean(xref_cbg[xref_cbg>0])
xtar_cbg_non0_mean = mean(xtar_cbg[xtar_cbg>0])
### get A & B
B = (log(xref_cpk_non0_mean) - log(xref_cbg_non0_mean)) / (log(xtar_cpk_non0_mean) - log(xtar_cbg_non0_mean))
A = xref_cbg_non0_mean / (xtar_cbg_non0_mean^B)
### s3norm tar sig
return(xtar^B * A)
}

### get highest frip
print('get highest frip')
dsig_linear_frip = apply(dsig_linear, 2, function(x) get_frip(x, 0.1))
ref_sig = dsig_linear[,dsig_linear_frip==max(dsig_linear_frip)]

dsig_scale = apply(dsig, 2, function(x) (x-mean(x))/sd(x)*median_sd + median_mean )

print('S3norm')
dsig_s3norm = apply(dsig_linear, 2, function(x) s3norm(ref_sig, x, 0.1) )
dsig_s3norm = log2(dsig_s3norm+0.001)

d_s3norm = cbind(d[,1:3], dsig_s3norm)
write.table(d_s3norm, 'rnaTPM.s3norm.txt', sep=' ', quote=F, col.names=T, row.names=F)
d_scale = cbind(d[,1:3], dsig_scale)
write.table(d_scale, 'rnaTPM.scale.txt', sep=' ', quote=F, col.names=T, row.names=F)




