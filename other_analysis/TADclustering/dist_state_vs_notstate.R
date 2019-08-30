### read read counts
d_rc = read.table('/storage/home/g/gzx103/group/HiC/mouse/processed/10kb/G1E-ER4all3.merged.chr19.10kb.matrix', header=F)

plot_region = c(1:dim(d_rc)[1])

d_rc_lim = d_rc
upperlim = 1000
lowerlim = 0
d_rc_lim[d_rc_lim>upperlim] = upperlim
d_rc_lim[d_rc_lim<lowerlim] = lowerlim

#colnames(d_rc_lim) = rep('', dim(d_rc_lim)[2])
#rownames(d_rc_lim) = rep('', dim(d_rc_lim)[1])


### get color matrix
d_rc_mat_ideas = d_rc
d_rc_mat_ideas[,] = 0


si=2
print(si)
bed = read.table(paste('mm10.10KB.chr19.s', si, '.bed', sep=''), header=F)
print(dim(bed))
d_rc_mat_ideas_s2 = d_rc_lim
d_rc_mat_ideas_s2[,] = 0
d_rc_mat_ideas_s2_binary = d_rc_mat_ideas_s2
plot_region_min = min(plot_region)
plot_region_max = max(plot_region)
plot_region_length = length(plot_region)
for (i in 1:dim(bed)[1]){
print(i)
bi = bed[i,4]
for (j in 1:dim(bed)[1]){
bj = bed[j,4]	
#if ((bi-plot_region_min+1>0) & (bj-plot_region_min+1>0) & (bi-plot_region_min+1<plot_region_length) & (bj-plot_region_min+1<plot_region_length)){
sig = d_rc_lim[bi,bj]
d_rc_mat_ideas_s2[bi,bj] = sig
d_rc_mat_ideas_s2_binary[bi,bj] = 1
#}
}
}



dist_mean_s2 = c()
dist_sd_s2 = c()
dist_n_s2 = c()
for (k in 1:200){
print(k)
sig_dist_vec_tmp = c()
for (j in k:(dim(d_rc_mat_ideas_s2_binary)[1])){
#j=1
d_rc_mat_ideas_s2_binary_ij = d_rc_mat_ideas_s2_binary[j-k+1,j]
if (d_rc_mat_ideas_s2_binary_ij==1){
sig_dist_vec_tmp = c(sig_dist_vec_tmp, d_rc_mat_ideas_s2[j-k+1,j])
}
}
dist_mean_s2[k] = mean(log2(sig_dist_vec_tmp+1))
dist_sd_s2[k] = sd(log2(sig_dist_vec_tmp+1))
dist_n_s2[k] = length(sig_dist_vec_tmp)
}
dist_mean_sd_s2 = cbind(dist_mean_s2, dist_sd_s2, dist_n_s2)
write.table(dist_mean_sd_s2, 'dist_mean_sd_s2.txt', quote=F, sep='\t', col.names=F, row.names=F)





d_rc_mat_ideas_snot2 = d_rc
d_rc_mat_ideas_snot2[,] = 0
d_rc_mat_ideas_snot2_binary = d_rc_mat_ideas_snot2
plot_region_min = min(plot_region)
plot_region_max = max(plot_region)
plot_region_length = length(plot_region)

si = 2
bed0 = read.table(paste('mm10.10KB.chr19.s', si, '.bed', sep=''), header=F)
bed_not = c()
for (si in c(1,3:14)){
print(si)
bed = read.table(paste('mm10.10KB.chr19.s', si, '.bed', sep=''), header=F)
print(dim(bed))
bed_not = rbind(bed_not, bed)
}


for (i in 1:dim(bed0)[1]){
print(i)
bi = bed0[i,4]
for (j in 1:dim(bed_not)[1]){
bj = bed_not[j,4]
sig = d_rc_lim[bi,bj]
d_rc_mat_ideas_snot2[bi,bj] = sig
d_rc_mat_ideas_snot2_binary[bi,bj] = 1
}
}



dist_mean_snot2 = c()
dist_sd_snot2 = c()
dist_n_snot2 = c()
#for (k in 1:(dim(d_rc_mat_ideas_snot2_binary)[1])){
for (k in 1:200){
print(k)
sig_dist_vec_tmp = c()
for (j in k:(dim(d_rc_mat_ideas_snot2_binary)[1])){
#j=1
d_rc_mat_ideas_s2_binary_ij = d_rc_mat_ideas_snot2_binary[j-k+1,j]
if (d_rc_mat_ideas_s2_binary_ij==1){
sig_dist_vec_tmp = c(sig_dist_vec_tmp, d_rc_mat_ideas_snot2[j-k+1,j])
}
}
dist_mean_snot2[k] = mean(log2(sig_dist_vec_tmp+1))
dist_sd_snot2[k] = sd(log2(sig_dist_vec_tmp+1))
dist_n_snot2[k] = length(sig_dist_vec_tmp)
}
dist_mean_sd_snot2 = cbind(dist_mean_snot2, dist_sd_snot2, dist_n_snot2)
#write.table(dist_mean_sd_snot2, 'dist_mean_sd_snot2.neg.txt', quote=F, sep='\t', col.names=F, row.names=F)




win_size = 200
pdf('signal_s2_vs_snot2.pdf', height=14)
par(mfrow=c(2,1))
plot(c(1:dim(dist_mean_sd_s2)[1])[1:win_size], dist_mean_sd_s2[,1][1:win_size], col='blue', log='', pch=16)
lines(c(1:dim(dist_mean_sd_s2)[1])[1:win_size], dist_mean_sd_s2[,1][1:win_size], col='blue')
points(c(1:dim(dist_mean_sd_s2)[1])[1:win_size], dist_mean_sd_snot2[,1][1:win_size], col='black', pch=16)
lines(c(1:dim(dist_mean_sd_s2)[1])[1:win_size], dist_mean_sd_snot2[,1][1:win_size], col='black')
dist_mean_sd_CP_s2 = cbind(dist_mean_sd_snot2[1:200,], dist_mean_sd_s2[1:200,])
t = apply(dist_mean_sd_CP_s2, 1, function(x) (x[4]-x[1])/sqrt(x[5]^2/x[6] + x[2]^2/x[3]) )  
df = apply(dist_mean_sd_CP_s2, 1, function(x) (x[5]/x[6]+x[2]/x[3])/(x[5]^4/(x[6]^2*(x[3]-1)) + x[2]^4/(x[3]^2*(x[6]-1))) )  
pt = apply(cbind(t, df), 1, function(x) 1*pt(x[1], x[2], lower=FALSE))
plot(c(1:dim(dist_mean_sd_s2)[1])[1:win_size], -log10(p.adjust(pt, 'fdr')[1:win_size]), col='blue', log='', pch=16, ylim = c(0, 16))
lines(c(1:dim(dist_mean_sd_s2)[1])[1:win_size], -log10(p.adjust(pt, 'fdr')[1:win_size]), col='blue')
abline(h=-log10(0.05))
dev.off()








