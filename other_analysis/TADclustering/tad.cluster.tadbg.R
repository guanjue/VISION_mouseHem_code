library(pheatmap)
library(mclust)
library(mixtools)
library(psych)


IDEAS = read.table('/storage/home/gzx103/group/projects/vision/pknorm_2_16lim_ref1mo_sample1mo_0424/get_eRP_tracks/ER4.state.bed', header=F)
head(IDEAS)
IDEAS_state_count_tad = table(IDEAS[,4])

### read info
d1 = read.table('tad.level.1.mat.p.txt', header=F)
d2 = read.table('tad.level.2.mat.p.txt', header=F)
d3 = read.table('tad.level.3.mat.p.txt', header=F)
d4 = read.table('tad.level.4.mat.p.txt', header=F)
d5 = read.table('tad.level.5.mat.p.txt', header=F)

d1c = read.table('tad.level.1.mat.c.txt', header=F)
d2c = read.table('tad.level.2.mat.c.txt', header=F)
d3c = read.table('tad.level.3.mat.c.txt', header=F)
d4c = read.table('tad.level.4.mat.c.txt', header=F)
d5c = read.table('tad.level.5.mat.c.txt', header=F)
d12345_info_c = rbind(d1c,d2c,d3c,d4c,d5c)
d12345_info_c0 = d12345_info_c[,-c(1:5)]
#IDEAS_state_count_tad = colSums(d12345_info_c0)


### merge TAD level
d12345_info0 = rbind(d1,d2,d3,d4,d5)
d12345_info1 = d12345_info0[order(d12345_info0[,2]),]
d12345_info = d12345_info1[order(d12345_info1[,1]),]

d12345 = d12345_info[,c(4,6:dim(d12345_info)[2])]
d12345_info_pk = d12345_info[,1:3]
d12345_info_pk_name = c()
for (i in 1:dim(d12345_info_pk)[1]){
	a = d12345_info_pk[i,]
	d12345_info_pk_name = rbind(d12345_info_pk_name, paste((a[1]), (a[2]), (a[3]), sep='_'))
}
d12345_info_pk_withname = cbind(d12345_info_pk, d12345_info_pk_name)
write.table(d12345_info_pk_withname, 'TAD_regions.txt', quote=FALSE, col.names=FALSE, row.names=FALSE, sep=' ')


### scale
d12345_s = d12345
d12345_s[,-1] = t(apply(d12345[,-1], 1, function(x) x/sum(x) ))
d12345_s[,-1] = t(apply(d12345[,-1], 1, function(x) x/IDEAS_state_count_tad*mean(IDEAS_state_count_tad) ))
d12345_s[,1] = ((d12345[,1]) - mean(d12345[,1]))/sd(d12345[,1]) * sd(as.matrix(d12345[,-1])) + mean(as.matrix(d12345[,-1]))

pdf('all_p_hist.pdf', height=200)
par(mfrow=c(31,1))
for (i in 1:28){
	hist(log10(d12345_s+0.1)[,i], breaks=50)
}
hist(as.matrix(log10(d12345_s[,-1]+0.1)), breaks=50)
hist(as.matrix(log10(d12345_s[,-1]+0.1)[log10(d12345_s[,-1]+0.1)>-1]), breaks=50)
hist(as.matrix((d12345_s[,-1]+0)), breaks=1000, xlim=c(0,5))
dev.off()

### get IDEAS input
log10_histp = (d12345_s[,-1]+0)

dir.create('IDEAS_input')
for (i in 0:26){
	log10_histp_i = log10_histp[,i+1]
	write.table(log10_histp_i, paste('IDEAS_input/IDEAS_state_p.', i, '.txt', sep=''), quote=FALSE, col.names=FALSE, row.names=FALSE)
}


### get IDEAS state p correlation matrix
histp_cor_mat = cor(log10_histp)
#cor_dist = as.dist(1-histp_cor_mat)
histp_cor_dist_mat = cor2dist(histp_cor_mat)
cor_dist = as.dist(histp_cor_dist_mat)

### get JI distance
hc = hclust(cor_dist, "complete")
###
dist_within = c()
nclustr = 27
for (i in c(1:nclustr)){
print(i)
cluster_id = cutree(hc, i)
dist_within_i = 0
for (j in c(1:i)){
dist_within_j = sum(histp_cor_dist_mat[cluster_id==j, cluster_id==j])
dist_within_i = dist_within_i + dist_within_j
}
dist_within[i] = dist_within_i
}
### plot JI dist curve
pdf('IDEASp_dist.curve.pdf')
par(mfrow=c(1,1))
plot(1:nclustr, dist_within, log='')
lines(1:nclustr, dist_within)
dev.off()

### ncluster
ncluster = 10
cluster_id_used = cutree(hc, ncluster)
corlo_lim_scale_max = 1
corlo_lim_scale_min = -1
breaksList_scale = seq(corlo_lim_scale_min, corlo_lim_scale_max, by = 0.001)
my_colorbar_scale=colorRampPalette(c('blue', 'white', 'red'))(n = length(breaksList_scale))
colnames(histp_cor_mat) = c(0:26)
rownames(histp_cor_mat) = c(0:26)
pdf('tad.ideasp.cor.pdf', width = 7, height = 7)
pheatmap(as.matrix(histp_cor_mat), clustering_distance_cols=cor_dist, clustering_distance_rows=cor_dist, color=my_colorbar_scale, breaks = breaksList_scale, annotation_names_col = FALSE)
dev.off()
cluster_colors = c('#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', '#46f0f0', '#f032e6', '#bcf60c', '#fabebe', '#008080', '#e6beff', '#9a6324', '#fffac8', '#800000', '#aaffc3', '#808000', '#ffd8b1', '#000075', '#808080', '#ffffff', '#000000')
cluster_colors = cluster_colors[1:ncluster]
cluster_mat = as.matrix(cbind(cluster_id_used, cluster_id_used))
colnames(cluster_mat) = c(20,20)
pdf('tad.ideasp.cor.cluster.pdf', width = 7, height = 7)
pheatmap(cluster_mat, clustering_distance_rows=cor_dist, cluster_cols=FALSE, color=cluster_colors, annotation_names_col = FALSE)
dev.off()

cluster_ids = unique(cluster_id_used[hc$order])

### merge IDEASps
d12345_IDEASps = d12345[,-1]
d12345_IDEASps_merged = c()
IDEAS_state_count_tad_merged = c()
for (i in 1:ncluster){
	if (sum(cluster_id_used==cluster_ids[i])>1){
	d12345_IDEASps_tmp = rowSums(d12345_IDEASps[,cluster_id_used==cluster_ids[i]])
	IDEAS_state_count_tad_tmp = sum(IDEAS_state_count_tad[cluster_id_used==cluster_ids[i]])
	} else{
	d12345_IDEASps_tmp = d12345_IDEASps[,cluster_id_used==cluster_ids[i]]
	IDEAS_state_count_tad_tmp = IDEAS_state_count_tad[cluster_id_used==cluster_ids[i]]
	}
	d12345_IDEASps_merged = cbind(d12345_IDEASps_merged, d12345_IDEASps_tmp)
	IDEAS_state_count_tad_merged = c(IDEAS_state_count_tad_merged, IDEAS_state_count_tad_tmp)
}
### scale merge IDEASps
d12345_IDEASps_merged_s = d12345_IDEASps_merged#t(apply(d12345_IDEASps_merged, 1, function(x) x/IDEAS_state_count_tad_merged*mean(IDEAS_state_count_tad_merged) ))

d12345_IDEASps_merged_s_log10 = (d12345_IDEASps_merged_s+0)

dir.create('IDEAS_input_merged')
for (i in 1:ncluster){
	log10_histp_i = d12345_IDEASps_merged_s_log10[,i]
	write.table(log10_histp_i, paste('IDEAS_input_merged/IDEAS_state_p.', i, '.txt', sep=''), quote=FALSE, col.names=FALSE, row.names=FALSE)
}









### fit GMM
log10_histp = as.vector(as.matrix(log10(d12345_s[,-1]+0.1)))

set.seed(2018)
gmm_histp = densityMclust(log10_histp, G=2, prior=priorControl(functionName="defaultPrior", shrinkage=0.1))
cluster_id = gmm_histp$classification
cluster_mean = gmm_histp$parameters$mean
print('2nd GMM cluster means: ')
print(cluster_mean_0hr)
#
rainbow_cp_0hr = rev(rainbow(length(cluster_mean_0hr)))
#
pdf('gmm_histp.pdf', width=7, height=7)
plot(gmm_histp, what = "density", data = log10_histp, breaks = 50)
for (i in c(1:length(cluster_mean))){
        print(i)
        x_input = seq(-2, 10, 0.05)#seq(-10, 140)
        cp_i_mean = gmm_histp$parameters$mean[i]
        cp_i_sd = (gmm_histp$parameters$variance$sigmasq[1])^0.5
        cp_i_pro = gmm_histp$parameters$pro[i]
        lines(x_input, cp_i_pro * dnorm(x_input, mean=cp_i_mean, sd=cp_i_sd), col=rainbow_cp_0hr[i])
}
dev.off()

### get TAD index
log10_histp_index = as.matrix(log10(d12345_s[,-1]+0.1))
log10_histp_index[as.matrix(log10(d12345_s[,-1]+0.1))>-1] = 1
log10_histp_index[as.matrix(log10(d12345_s[,-1]+0.1))==-1] = 0
log10_histp_IS = apply(log10_histp_index, 1, function(x) paste(x, collapse='_'))

### get TAD thresh
index_count = table(log10_histp_IS)
pdf('IS_count.pdf')
hist(log2(index_count), breaks=50)
dev.off()
index_count_thresh = 3






### add colnames
colnames(d12345_s) = c('TAD', 0:26)
colnames(d12345) = c('TAD', 0:26)

### log10
d12345_s_log10 = log10(d12345_s+0.1)
d12345_s_s = d12345_s_log10[,-1]

### splot TAD levels
d12345_s_s1 = d12345_s_s[(1):(dim(d1)[1]),]
d12345_s_s2 = d12345_s_s[(1+dim(d1)[1]):(dim(d1)[1]+dim(d2)[1]),]
d12345_s_s3 = d12345_s_s[(1+dim(d1)[1]+dim(d2)[1]):(dim(d1)[1]+dim(d2)[1]+dim(d3)[1]),]
d12345_s_s4 = d12345_s_s[(1+dim(d1)[1]+dim(d2)[1]+dim(d3)[1]):(dim(d1)[1]+dim(d2)[1]+dim(d3)[1]+dim(d4)[1]),]
d12345_s_s5 = d12345_s_s[(1+dim(d1)[1]+dim(d2)[1]+dim(d3)[1]+dim(d4)[1]):(dim(d1)[1]+dim(d2)[1]+dim(d3)[1]+dim(d4)[1]+dim(d5)[1]),]

### Mclust
set.seed(2018)
fit1 = Mclust(d12345_s_s1, G=5)
fit2 = Mclust(d12345_s_s2, G=5)
fit3 = Mclust(d12345_s_s3, G=5)
fit4 = Mclust(d12345_s_s4, G=5)
fit5 = Mclust(d12345_s_s5, G=5)


### sort 
d12345_s_s = c()
d12345_s_s = rbind(d12345_s_s, d12345_s_s1[order(fit1$classification),])
d12345_s_s = rbind(d12345_s_s, d12345_s_s2[order(fit2$classification),])
d12345_s_s = rbind(d12345_s_s, d12345_s_s3[order(fit3$classification),])
d12345_s_s = rbind(d12345_s_s, d12345_s_s4[order(fit4$classification),])
d12345_s_s = rbind(d12345_s_s, d12345_s_s5[order(fit5$classification),])

### plot
corlo_lim_scale_max = max(d12345_s_s)
corlo_lim_scale_min = min(d12345_s_s)
breaksList_scale = seq(corlo_lim_scale_min, corlo_lim_scale_max, by = 0.001)
my_colorbar_scale=colorRampPalette(c('white', 'red'))(n = length(breaksList_scale))

pdf('tad.cluster.scale.log10.5levels.pdf', width = 7, height = 7)
pheatmap(as.matrix(d12345_s_s), color=my_colorbar_scale, breaks = breaksList_scale, annotation_names_col = FALSE, cluster_cols=T, cluster_rows=F)
dev.off()


set.seed(2018)
fit0 = Mclust(d12345_s_s, G=50)
d12345_s_s = d12345_s_s[order(fit0$classification),]

### plot
corlo_lim_scale_max = max(d12345_s_s)
corlo_lim_scale_min = min(d12345_s_s)
breaksList_scale = seq(corlo_lim_scale_min, corlo_lim_scale_max, by = 0.001)
my_colorbar_scale=colorRampPalette(c('white', 'red'))(n = length(breaksList_scale))

pdf('tad.cluster.scale.log10.1levels.pdf', width = 7, height = 7)
pheatmap(as.matrix(d12345_s_s), color=my_colorbar_scale, breaks = breaksList_scale, annotation_names_col = FALSE, cluster_cols=T, cluster_rows=F)
dev.off()



d12345_s_s = d12345_s[,-1]
corlo_lim_scale_max = 5
corlo_lim_scale_min = min(d12345_s_s)
breaksList_scale = seq(corlo_lim_scale_min, corlo_lim_scale_max, by = 0.001)
my_colorbar_scale=colorRampPalette(c('white', 'red'))(n = length(breaksList_scale))

pdf('tad.cluster.scale.pdf', width = 7, height = 7)
pheatmap(as.matrix(d12345_s_s), color=my_colorbar_scale, breaks = breaksList_scale, annotation_names_col = FALSE, cluster_cols=T, cluster_rows=T)
dev.off()


corlo_lim_scale_max = 10
corlo_lim_scale_min = min(scale(d12345[,-1]))
breaksList_scale = seq(corlo_lim_scale_min, corlo_lim_scale_max, by = 0.001)
my_colorbar_scale=colorRampPalette(c('white', 'red'))(n = length(breaksList_scale))
pdf('tad.cluster.pdf', width = 7, height = 7)
pheatmap(as.matrix(scale(d12345[,-1])), color=my_colorbar_scale, breaks = breaksList_scale, annotation_names_col = FALSE, cluster_cols=T, cluster_rows=T)
dev.off()










