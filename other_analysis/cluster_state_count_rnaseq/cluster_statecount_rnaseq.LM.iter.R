library(randomForest)
library(pheatmap)
library(e1071)
library(neuralnet)
library(LSD)

########### R2
R2 = function(d_IP_ts, qPCR_expand){
        r2 = 1-sum((qPCR_expand-d_IP_ts)^2) / sum((qPCR_expand-mean(qPCR_expand))^2)
        return(r2)
}
###########

sigmat = read.table('~/group/projects/vision/pknorm_2_16lim_ref1mo_sample1mo_0424/pknorm_2_16lim_ref1mo_0424_lesshet.para0')
sigmat1 = as.matrix(sigmat[,2:9]/sigmat[,1])
print(head(sigmat1))
fit=hclust(dist(sigmat1),method="ward.D2")


### get gene group
get_gene_group = function(gene_exp){
	### get mean and sd
	m = apply(gene_exp, 1, mean)
	s = apply(gene_exp, 1, sd)
	tt3 = which(m>-4 & s>2)
	tt4 = which(m>-4 & s<=2)
	tt2 = which(m<=-4 & s>2)
	tt1 = which(m<=-4 & s<=2)
	### get gene_group_label
	label_all_genegroup1 = rep(0, length(m))
	label_all_genegroup1[tt3] = 3
	label_all_genegroup1[tt4] = 4
	label_all_genegroup1[tt2] = 2
	label_all_genegroup1[tt1] = 1
	return(label_all_genegroup1)
}

### get predictor (state count) colnames
colnames_state = c()
for (i in 0:26){
	colnames_state = c(colnames_state, paste('D',i,sep=''))
}
for (i in 0:26){
	colnames_state = c(colnames_state, paste('P',i,sep=''))
}
colnames_state = c(colnames_state, 'Gene_Group')



cell_num = c(1:12)
chr_num_all = c(1:19, 'X')
cell_num = c(1:3)
chr_num_all = c(1:3, 'X')

cell_num_test = 1
chr_num_test = '1'

LM_r2 = c()
LM_P_r2 = c()
LM_D_r2 = c()
GG_r2_benchmark = c()
GG_r2 = c()
GG_P_r2 = c()
GG_D_r2 = c()

for (tcell in cell_num){
for (tchr in chr_num_all){

cell_num_test = tcell
chr_num_test = tchr

gene_exp = c()
gene_exp_genegroup = c()
dr0_all_mat = c()
gene_exp0 = c()
gene_exp_genegroup0 = c()
dr0_all_mat0 = c()

### get input matrix
for (i in c(1:4)){
gene_exp_tmp = c()
dr0_all_tmp = c()
rt_y = c()
rt_x = c()
rt_x0 = c()
### get dif training chr data
for (k in chr_num_all[chr_num_all!=chr_num_test]){
	load(paste('../vision_rna_tss2k_ccreunit.chr', k, '.', i, '.', 1, '.Rdat', sep=''))
	rt_y = c(rt_y, rt$y)
	rt_x = rbind(rt_x, rt$x)
	rt_x0 = rbind(rt_x0, rt$x0)
}
### get cell type mat
for (j in cell_num[cell_num!=cell_num_test]){
	used_id_ct = seq(1,length(rt_y),by=12) + (j-1)
	gene_exp_tmp = cbind(gene_exp_tmp, rt_y[used_id_ct])
	dr0_all_tmp = cbind(dr0_all_tmp, cbind(rt_x[used_id_ct,], rt_x0[used_id_ct,]) )
}
### get reordered gene group id
gene_exp_genegroup = c(gene_exp_genegroup, get_gene_group(gene_exp_tmp))
gene_exp = rbind(gene_exp, gene_exp_tmp)
dr0_all_mat = rbind(dr0_all_mat, dr0_all_tmp)
gene_exp0 = c(gene_exp0, rt_y)
dr0_all_mat0 = rbind(dr0_all_mat0,  cbind(rt_x, rt_x0))
}


### get training input matrix
gene_exp_vec = c()
gene_exp_genegroup_vec = c()
dr0_all = c()
for (i in c(1:dim(gene_exp)[2])){
	### vector gene exp mat
	gene_exp_vec = c(gene_exp_vec, gene_exp[,i])
	gene_exp_genegroup_vec = c(gene_exp_genegroup_vec, gene_exp_genegroup)
	### split DP col
	used_col = c(1:54) + (i-1) * 54
	dr0_all = rbind(dr0_all, dr0_all_mat[,used_col])
}
dr0_all = as.data.frame(cbind(dr0_all, gene_exp_vec))
colnames(dr0_all) = colnames_state


### get testing cell chr gene group
gene_exp_genegroup_test = c()
gene_exp_genegroup_forbenchmark = c()
for (i in c(1:4)){
load(paste('../vision_rna_tss2k_ccreunit.chr', chr_num_test, '.', i, '.', cell_num_test, '.Rdat', sep=''))
gene_exp_tmp = c()
for (j in cell_num[cell_num!=cell_num_test]){
	used_id_ct = seq(1,length(rt$y),12) + (j-1)
	gene_exp_tmp = cbind(gene_exp_tmp, rt$y[used_id_ct])
}
gene_exp_genegroup_test = c(gene_exp_genegroup_test, get_gene_group(gene_exp_tmp))
gene_exp_genegroup_forbenchmark = rbind(gene_exp_genegroup_forbenchmark, gene_exp_tmp)
}

### get predictor mat and exp vec of testing data
gene_exp_test = c()
dr0_all_mat_test = c()
### get input matrix
for (i in c(1:4)){
load(paste('../vision_rna_tss2k_ccreunit.chr', chr_num_test, '.', i, '.', cell_num_test, '.Rdat', sep=''))
gene_exp_tmp = c()
dr0_all_tmp = c()
for (j in c(cell_num_test)){
	used_id_ct = seq(1,length(rt$y),12) + (j-1)
	gene_exp_tmp = cbind(gene_exp_tmp, rt$y[used_id_ct])
	dr0_all_tmp = cbind(dr0_all_tmp, cbind(rt$x[used_id_ct,], rt$x0[used_id_ct,]) )
}
gene_exp_test = rbind(gene_exp_test, gene_exp_tmp)
dr0_all_mat_test = rbind(dr0_all_mat_test,  dr0_all_tmp)
}

### get testing input
gene_exp_vec_test = gene_exp_test
dr0_all_testing = dr0_all_mat_test
dr0_all_testing = as.data.frame(cbind(dr0_all_testing, gene_exp_vec_test))
colnames(dr0_all_testing) = colnames_state


############### 
### prediction
### All
fit = lm(Gene_Group ~ .-1, data=dr0_all)
fit_R2 = summary(fit)$r.squared
fit_R2
exp_fitpre = predict(fit, newdata = dr0_all_testing[,1:54])
fit_R2pre = R2(exp_fitpre, dr0_all_testing[,55])
fit_R2pre
LM_r2 = c(LM_r2, fit_R2pre)
0.4996497
### Proximal
fit_P = lm(Gene_Group ~ .-1, data=dr0_all[,c(28:55)])
fit_R2_P = summary(fit_P)$r.squared
fit_R2_P
exp_fitpre_P = predict(fit_P, newdata = dr0_all_testing[,c(28:54)])
fit_R2pre_P = R2(exp_fitpre_P, dr0_all_testing[,55])
fit_R2pre_P
LM_P_r2 = c(LM_P_r2, fit_R2pre_P)
0.3540576
### Distal
fit_D = lm(Gene_Group ~ .-1, data=dr0_all[,c(1:27,55)])
fit_R2_D = summary(fit_D)$r.squared
fit_R2_D
exp_fitpre_D = predict(fit_D, newdata = dr0_all_testing[,c(1:27)])
fit_R2pre_D = R2(exp_fitpre_D, dr0_all_testing[,55])
fit_R2pre_D
LM_D_r2 = c(LM_D_r2, fit_R2pre_D)
0.2982498

png(paste('test.LM.chr', chr_num_test, '.', cell_num_test, '.png', sep=''), width=1500, height=500)
par(mfrow=c(1,3)) 
plot(exp_fitpre, dr0_all_testing[,55], xlim = c(-10, 15), ylim = c(-10, 15), main=round(R2(exp_fitpre, dr0_all_testing[,55]), 3))
abline(0,1, col='red')
plot(exp_fitpre_P, dr0_all_testing[,55], xlim = c(-10, 15), ylim = c(-10, 15), main=round(R2(exp_fitpre_P, dr0_all_testing[,55]), 3))
abline(0,1, col='red')
plot(exp_fitpre_D, dr0_all_testing[,55], xlim = c(-10, 15), ylim = c(-10, 15), main=round(R2(exp_fitpre_D, dr0_all_testing[,55]), 3))
abline(0,1, col='red')
dev.off()



### All GG
exp_fitpre_gg_all = c()
exp_fitpre_gg_all_benchmark = c()
dr0_all_testing_gg_all = c()
gene_exp_genegroup_vec_testing_new = c()

for (gg in c(1:4)){
print(gg)
print(sum(gene_exp_genegroup_vec==gg))
### get gene group training
dr0_all_gg = dr0_all[gene_exp_genegroup_vec==gg,]
### fit LM
fit_gg = lm(Gene_Group ~ .-1, data=dr0_all_gg)
fit_R2_gg = summary(fit_gg)$r.squared
print(fit_R2_gg)
### get gene group testing
dr0_all_testing_gg = dr0_all_testing[gene_exp_genegroup_test==gg,]
### bench mark
exp_fitpre_gg_bm_mean = mean(gene_exp_genegroup_forbenchmark[gene_exp_genegroup_test==gg,])
exp_fitpre_gg_bm = rep(exp_fitpre_gg_bm_mean, sum(gene_exp_genegroup_test==gg))
exp_fitpre_gg_all_benchmark = c(exp_fitpre_gg_all_benchmark, exp_fitpre_gg_bm)
### predict 
exp_fitpre_gg = predict(fit_gg, dr0_all_testing_gg)
### get R2
fit_R2pre_gg = R2(exp_fitpre_gg, dr0_all_testing_gg[,55])
### get pred vec
exp_fitpre_gg_all = c(exp_fitpre_gg_all, exp_fitpre_gg)
### get obs vec
dr0_all_testing_gg_all = c(dr0_all_testing_gg_all, dr0_all_testing_gg[,55])
### get gene group vec
gene_exp_genegroup_vec_testing_new = c(gene_exp_genegroup_vec_testing_new, gene_exp_genegroup_test[gene_exp_genegroup_test==gg])
print(fit_R2pre_gg)
}
fitpre_LM_R2_gg = R2(exp_fitpre_gg_all, dr0_all_testing_gg_all)
fitpre_LM_R2_gg
fitpre_LM_R2_gg_benchmark = R2(exp_fitpre_gg_all_benchmark, dr0_all_testing_gg_all)
fitpre_LM_R2_gg_benchmark
GG_r2_benchmark = c(GG_r2_benchmark, fitpre_LM_R2_gg_benchmark)
GG_r2 = c(GG_r2, fitpre_LM_R2_gg)
0.7456658
0.7055312
0.7024672

### All GG _P
exp_fitpre_gg_all_P = c()
dr0_all_testing_gg_all_P = c()
gene_exp_genegroup_vec_testing_new_P = c()

for (gg in c(1:4)){
print(gg)
print(sum(gene_exp_genegroup_vec==gg))
### get gene group training
dr0_all_gg = dr0_all[gene_exp_genegroup_vec==gg,]
### fit LM
fit_gg = lm(Gene_Group ~ .-1, data=dr0_all_gg[,28:55])
fit_R2_gg = summary(fit_gg)$r.squared
print(fit_R2_gg)
### get gene group testing
dr0_all_testing_gg = dr0_all_testing[gene_exp_genegroup_test==gg,]
### predict 
exp_fitpre_gg = predict(fit_gg, dr0_all_testing_gg[,28:54])
### get R2
fit_R2pre_gg = R2(exp_fitpre_gg, dr0_all_testing_gg[,55])
### get pred vec
exp_fitpre_gg_all_P = c(exp_fitpre_gg_all_P, exp_fitpre_gg)
### get obs vec
dr0_all_testing_gg_all_P = c(dr0_all_testing_gg_all_P, dr0_all_testing_gg[,55])
### get gene group vec
gene_exp_genegroup_vec_testing_new_P = c(gene_exp_genegroup_vec_testing_new, gene_exp_genegroup_test[gene_exp_genegroup_test==gg])
print(fit_R2pre_gg)
}
fitpre_LM_R2_gg_P = R2(exp_fitpre_gg_all_P, dr0_all_testing_gg_all_P)
fitpre_LM_R2_gg_P
GG_P_r2 = c(GG_P_r2, fitpre_LM_R2_gg_P)
0.6869577

### All GG _D
exp_fitpre_gg_all_D = c()
dr0_all_testing_gg_all_D = c()
gene_exp_genegroup_vec_testing_new_D = c()

for (gg in c(1:4)){
print(gg)
print(sum(gene_exp_genegroup_vec==gg))
### get gene group training
dr0_all_gg = dr0_all[gene_exp_genegroup_vec==gg,]
### fit LM
fit_gg = lm(Gene_Group ~ .-1, data=dr0_all_gg[,c(1:27,55)])
fit_R2_gg = summary(fit_gg)$r.squared
print(fit_R2_gg)
### get gene group testing
dr0_all_testing_gg = dr0_all_testing[gene_exp_genegroup_test==gg,]
### predict 
exp_fitpre_gg = predict(fit_gg, dr0_all_testing_gg[,c(1:27)])
### get R2
fit_R2pre_gg = R2(exp_fitpre_gg, dr0_all_testing_gg[,55])
### get pred vec
exp_fitpre_gg_all_D = c(exp_fitpre_gg_all_D, exp_fitpre_gg)
### get obs vec
dr0_all_testing_gg_all_D = c(dr0_all_testing_gg_all_D, dr0_all_testing_gg[,55])
### get gene group vec
gene_exp_genegroup_vec_testing_new_D = c(gene_exp_genegroup_vec_testing_new, gene_exp_genegroup_test[gene_exp_genegroup_test==gg])
print(fit_R2pre_gg)
}
fitpre_LM_R2_gg_D = R2(exp_fitpre_gg_all_D, dr0_all_testing_gg_all_D)
fitpre_LM_R2_gg_D
GG_D_r2 = c(GG_D_r2, fitpre_LM_R2_gg_D)
0.7239277


### plot GG
png(paste('test.GG.chr', chr_num_test, '.', cell_num_test, '.png', sep=''), width=1500, height=500)
par(mfrow=c(1,3)) 
plot(exp_fitpre_gg_all, dr0_all_testing_gg_all, xlim = c(-10, 15), ylim = c(-10, 15), main=round(R2(exp_fitpre_gg_all, dr0_all_testing_gg_all),3))
abline(0,1, col='red')
plot(exp_fitpre_gg_all_P, dr0_all_testing_gg_all_P, xlim = c(-10, 15), ylim = c(-10, 15), main=round(R2(exp_fitpre_gg_all_P, dr0_all_testing_gg_all_P),3))
abline(0,1, col='red')
plot(exp_fitpre_gg_all_D, dr0_all_testing_gg_all_D, xlim = c(-10, 15), ylim = c(-10, 15), main=round(R2(exp_fitpre_gg_all_D, dr0_all_testing_gg_all_D),3))
abline(0,1, col='red')
dev.off()

### plot split
png(paste('test.GG.split.chr', chr_num_test, '.', cell_num_test, '.png', sep=''), width=2000, height=500)
par(mfrow=c(1,4)) 
for (i in c(1:4)){
used_id_gg_plot = gene_exp_genegroup_vec_testing_new==i
plot(exp_fitpre_gg_all[used_id_gg_plot], dr0_all_testing_gg_all[used_id_gg_plot], xlim = c(-10, 15), ylim = c(-10, 15), main=round(R2(exp_fitpre_gg_all[used_id_gg_plot], dr0_all_testing_gg_all[used_id_gg_plot]),3))
abline(0,1, col='red')
}
dev.off()

png(paste('test.GG.split.P.chr', chr_num_test, '.', cell_num_test, '.png', sep=''), width=2000, height=500)
par(mfrow=c(1,4)) 
for (i in c(1:4)){
used_id_gg_plot = gene_exp_genegroup_vec_testing_new==i
plot(exp_fitpre_gg_all_P[used_id_gg_plot], dr0_all_testing_gg_all[used_id_gg_plot], xlim = c(-10, 15), ylim = c(-10, 15), main=round(R2(exp_fitpre_gg_all_P[used_id_gg_plot], dr0_all_testing_gg_all[used_id_gg_plot]),3))
abline(0,1, col='red')
}
dev.off()
png(paste('test.GG.split.D.chr', chr_num_test, '.', cell_num_test, '.png', sep=''), width=2000, height=500)
par(mfrow=c(1,4)) 
for (i in c(1:4)){
used_id_gg_plot = gene_exp_genegroup_vec_testing_new==i
plot(exp_fitpre_gg_all_D[used_id_gg_plot], dr0_all_testing_gg_all[used_id_gg_plot], xlim = c(-10, 15), ylim = c(-10, 15), main=round(R2(exp_fitpre_gg_all_D[used_id_gg_plot], dr0_all_testing_gg_all[used_id_gg_plot]),3))
abline(0,1, col='red')
}
dev.off()

print(paste('chr', chr_num_test, '.cell', cell_num_test, '.finished', sep=''))

}
}

write.table(LM_r2, 'LM_r2.txt', quote=F, col.names=F, row.names=F, sep='\t')
write.table(LM_P_r2, 'LM_P_r2.txt', quote=F, col.names=F, row.names=F, sep='\t')
write.table(LM_D_r2, 'LM_D_r2.txt', quote=F, col.names=F, row.names=F, sep='\t')
write.table(GG_r2_benchmark, 'GG_r2_benchmark.txt', quote=F, col.names=F, row.names=F, sep='\t')
write.table(GG_r2, 'GG_r2.txt', quote=F, col.names=F, row.names=F, sep='\t')
write.table(GG_P_r2, 'GG_P_r2.txt', quote=F, col.names=F, row.names=F, sep='\t')
write.table(GG_D_r2, 'GG_D_r2.txt', quote=F, col.names=F, row.names=F, sep='\t')

png('r2_box.png')
R2_mat = as.matrix(cbind(rep(0, length(LM_P_r2)), LM_P_r2, LM_D_r2, LM_r2, GG_r2_benchmark, GG_P_r2, GG_D_r2, GG_r2))
boxplot(R2_mat)
for (i in c(1:dim(R2_mat)[1])){
	#lines(c(1:(dim(cor_mat_t)[2])), cor_mat_t[i,], col=rgb(255/255,0/255,0/255,alpha=0.5) )
	points(c(1:(dim(R2_mat)[2])), R2_mat[i,], col='black', pch=20 )
}
dev.off()



### plot mean vs sd
m = apply(gene_exp, 1, mean)
s = apply(gene_exp, 1, sd)
library(mclust)

m_vs_sd = (cbind(m, s))
BIC <- mclustBIC(m_vs_sd)
summary(BIC)
mod1 <- Mclust(scale(m_vs_sd), G = 4)
summary(mod1, parameters = TRUE)
c = mod1$classification

png('test.m_vs_sd.png', width=1000, height=1000)
par(mfrow=c(2,2))
plot(m, s)
abline(v=-4, col='red')
abline(h=2, col='red')
plot(m, s)
points(m[c==1], s[c==1], col='red')
points(m[c==2], s[c==2],col='blue')
points(m[c==3], s[c==3],col='green')
points(m[c==4], s[c==4],col='gray')
abline(v=-4, col='red')
abline(h=2, col='red')
hist(m, breaks=50)
abline(v=-4, col='red')
box()
hist(s, breaks=50)
abline(v=2, col='red')
box()
dev.off()

