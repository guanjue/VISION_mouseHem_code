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


gene_exp = c()
gene_exp_genegroup = c()
dr0_all_mat = c()
gene_exp0 = c()
gene_exp_genegroup0 = c()
dr0_all_mat0 = c()
cell_num = c(1:12)
chr_num_all = c(1:19, 'X')
cell_num_test = 2
chr_num_test = '1'
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


### get training input
gene_exp_vec = c()
gene_exp_genegroup_vec = c()
dr0_all = c()
for (i in c(1:dim(gene_exp)[2])){
gene_exp_vec = c(gene_exp_vec, gene_exp[,i])
gene_exp_genegroup_vec = c(gene_exp_genegroup_vec, gene_exp_genegroup)
used_col = c(1:54) + (i-1) * 54
dr0_all = rbind(dr0_all, dr0_all_mat[,used_col])
}
dr0_all = as.data.frame(cbind(dr0_all, gene_exp_vec))
colnames(dr0_all) = colnames_state



gene_exp_genegroup_test = c()
gene_exp0_test = c()
dr0_all_mat0_test = c()
test_chr = 1
test_cell = c(1:12)[-trainingcol]
#test_cell = c(1:11)
### get input matrix
for (i in c(1:4)){
load(paste('../vision_rna_tss2k_ccreunit.chr', chr_num_test, '.', i, '.', 1, '.Rdat', sep=''))
gene_exp_tmp = c()
for (j in trainingcol){
used_id_ct = seq(1,length(rt$y),12) + (j-1)
gene_exp_tmp = cbind(gene_exp_tmp, rt$y[used_id_ct])
}
gene_exp_genegroup_test = c(gene_exp_genegroup_test, get_gene_group(gene_exp_tmp))
gene_exp0_test = c(gene_exp0_test, rt$y)
dr0_all_mat0_test = rbind(dr0_all_mat0_test, cbind(rt$x, rt$x0))
}

gene_exp_test = c()
dr0_all_mat_test = c()
### get input matrix
for (i in c(1:4)){
load(paste('../vision_rna_tss2k_ccreunit.chr', chr_num_test, '.', i, '.', 1, '.Rdat', sep=''))
gene_exp_tmp = c()
dr0_all_tmp = c()
for (j in c(test_cell)){
used_id_ct = seq(1,length(rt$y),12) + (j-1)
gene_exp_tmp = cbind(gene_exp_tmp, rt$y[used_id_ct])
dr0_all_tmp = cbind(dr0_all_tmp, cbind(rt$x[used_id_ct,], rt$x0[used_id_ct,]) )
}
gene_exp_test = rbind(gene_exp_test, gene_exp_tmp)
dr0_all_mat_test = rbind(dr0_all_mat_test,  dr0_all_tmp)
}

### get testing input
gene_exp_vec_test = c()
dr0_all_testing = c()
for (i in c(1:dim(gene_exp_test)[2])){
gene_exp_vec_test = c(gene_exp_vec_test, gene_exp_test[,i])
used_col = c(1:54) + (i-1) * 54
dr0_all_testing = rbind(dr0_all_testing, dr0_all_mat_test[,used_col])
}
dr0_all_testing = as.data.frame(cbind(dr0_all_testing, gene_exp_vec_test))
colnames(dr0_all_testing) = colnames_state





### All
fit111 = lm(Gene_Group ~ ., data=dr0_all)
fit_R2111 = summary(fit111)$r.squared
fit_R2111
fit = lm(Gene_Group ~ .-1, data=dr0_all)
fit_R2 = summary(fit)$r.squared
fit_R2
exp_fitpre = predict(fit, newdata = dr0_all_testing[,1:54])
fit_R2pre = R2(exp_fitpre, dr0_all_test[,55])
fit_R2pre
0.4701275
### Proximal
fit_P = lm(Gene_Group ~ .-1, data=dr0_all[,c(28:55)])
fit_R2_P = summary(fit_P)$r.squared
fit_R2_P
exp_fitpre_P = predict(fit_P, newdata = dr0_all_testing[,c(28:54)])
fit_R2pre_P = R2(exp_fitpre_P, gene_exp_test)
fit_R2pre_P
0.3971424
### Distal
fit_D = lm(Gene_Group ~ .-1, data=dr0_all[,c(1:27,55)])
fit_R2_D = summary(fit_D)$r.squared
fit_R2_D
exp_fitpre_D = predict(fit_D, newdata = dr0_all_testing[,c(1:27)])
fit_R2pre_D = R2(exp_fitpre_D, gene_exp_test)
fit_R2pre_D
0.1202532

png('test.LM.png', width=1500, height=500)
par(mfrow=c(1,3)) 
plot(exp_fitpre, gene_exp_test, xlim = c(-10, 15), ylim = c(-10, 15), main=round(R2(exp_fitpre, gene_exp_test), 3))
abline(0,1, col='red')
plot(exp_fitpre_P, gene_exp_test, xlim = c(-10, 15), ylim = c(-10, 15), main=round(R2(exp_fitpre_P, gene_exp_test), 3))
abline(0,1, col='red')
plot(exp_fitpre_D, gene_exp_test, xlim = c(-10, 15), ylim = c(-10, 15), main=round(R2(exp_fitpre_D, gene_exp_test), 3))
abline(0,1, col='red')
dev.off()



### All GG
exp_fitpre_gg_all = c()
dr0_all_testing_gg_all = c()
exp_fitpre_gg_all_training = c()
dr0_all_testing_gg_all_training = c()
gene_exp_genegroup_vec_new = c()
gene_exp_genegroup_vec_testing_new = c()
dr0_all_new = c()
for (gg in c(1:4)){
	print(gg)
	print(sum(gene_exp_genegroup_vec==gg))
dr0_all_gg = dr0_all[gene_exp_genegroup_vec==gg,]
dr0_all_new = rbind(dr0_all_new, dr0_all_gg)
fit_gg = lm(Gene_Group ~ .-1, data=dr0_all_gg[,])
exp_fitpre_gg_all_training = c(exp_fitpre_gg_all_training, fit_gg$fitted.values)
dr0_all_testing_gg_all_training = c(dr0_all_testing_gg_all_training, dr0_all_gg[,55])
fit_R2_gg = summary(fit_gg)$r.squared
print(fit_R2_gg)
dr0_all_testing_gg = dr0_all_testing[gene_exp_genegroup_test==gg,]
exp_fitpre_gg = predict(fit_gg, dr0_all_testing_gg[,])
fit_R2pre_gg = R2(exp_fitpre_gg, dr0_all_testing_gg[,55])
exp_fitpre_gg_all = c(exp_fitpre_gg_all, exp_fitpre_gg)
dr0_all_testing_gg_all = c(dr0_all_testing_gg_all, dr0_all_testing_gg[,55])
gene_exp_genegroup_vec_new = c(gene_exp_genegroup_vec_new, gene_exp_genegroup_vec[gene_exp_genegroup_vec==gg])
gene_exp_genegroup_vec_testing_new = c(gene_exp_genegroup_vec_testing_new, gene_exp_genegroup_test[gene_exp_genegroup_test==gg])
print(fit_R2pre_gg)
}
fitpre_LM_R2_gg = R2(exp_fitpre_gg_all, dr0_all_testing_gg_all)
fitpre_LM_R2_gg
0.7571352
### All GG _P
exp_fitpre_gg_all_P = c()
dr0_all_testing_gg_all_P = c()
for (gg in c(1:4)){
	print(gg)
	print(sum(gene_exp_genegroup_vec==gg))
dr0_all_gg = dr0_all[gene_exp_genegroup_vec==gg,]
fit_gg_P = lm(Gene_Group ~ .-1, data=dr0_all_gg[,28:55])
fit_R2_gg_P = summary(fit_gg_P)$r.squared
print(fit_R2_gg_P)
dr0_all_testing_gg_P = dr0_all_testing[gene_exp_genegroup_test==gg,]
exp_fitpre_gg_P = predict(fit_gg_P, dr0_all_testing_gg_P[,28:55])
fit_R2pre_gg_P = R2(exp_fitpre_gg_P, dr0_all_testing_gg_P[,55])
exp_fitpre_gg_all_P = c(exp_fitpre_gg_all_P, exp_fitpre_gg_P)
dr0_all_testing_gg_all_P = c(dr0_all_testing_gg_all_P, dr0_all_testing_gg_P[,55])
print(fit_R2pre_gg_P)
}
fitpre_LM_R2_gg_P = R2(exp_fitpre_gg_all_P, dr0_all_testing_gg_all_P)
fitpre_LM_R2_gg_P
0.7169234
### All GG _D
exp_fitpre_gg_all_D = c()
dr0_all_testing_gg_all_D = c()
for (gg in c(1:4)){
	print(gg)
	print(sum(gene_exp_genegroup_vec==gg))
dr0_all_gg = dr0_all[gene_exp_genegroup_vec==gg,]
fit_gg_D = lm(Gene_Group ~ .-1, data=dr0_all_gg[,c(1:27,55)])
fit_R2_gg_D = summary(fit_gg_D)$r.squared
print(fit_R2_gg_D)
dr0_all_testing_gg_D = dr0_all_testing[gene_exp_genegroup_test==gg,]
exp_fitpre_gg_D = predict(fit_gg_D, dr0_all_testing_gg_D[,c(1:27,55)])
fit_R2pre_gg_D = R2(exp_fitpre_gg_D, dr0_all_testing_gg_D[,55])
exp_fitpre_gg_all_D = c(exp_fitpre_gg_all_D, exp_fitpre_gg_D)
dr0_all_testing_gg_all_D = c(dr0_all_testing_gg_all_D, dr0_all_testing_gg_D[,55])
print(fit_R2pre_gg_D)
}
fitpre_LM_R2_gg_D = R2(exp_fitpre_gg_all_D, dr0_all_testing_gg_all_D)
fitpre_LM_R2_gg_D
0.7450707

### plot GG
png('test.GG.png', width=1500, height=500)
par(mfrow=c(1,3)) 
plot(exp_fitpre_gg_all, dr0_all_testing_gg_all, xlim = c(-10, 15), ylim = c(-10, 15), main=round(R2(exp_fitpre_gg_all, dr0_all_testing_gg_all),3))
abline(0,1, col='red')
plot(exp_fitpre_gg_all_P, dr0_all_testing_gg_all_P, xlim = c(-10, 15), ylim = c(-10, 15), main=round(R2(exp_fitpre_gg_all_P, dr0_all_testing_gg_all_P),3))
abline(0,1, col='red')
plot(exp_fitpre_gg_all_D, dr0_all_testing_gg_all_D, xlim = c(-10, 15), ylim = c(-10, 15), main=round(R2(exp_fitpre_gg_all_D, dr0_all_testing_gg_all_D),3))
abline(0,1, col='red')
dev.off()

### plot split
png('test.GG.split.LM.png', width=2000, height=500)
par(mfrow=c(1,4)) 
for (i in c(1:4)){
used_id_gg_plot = gene_exp_genegroup_vec_testing_new==i
plot(exp_fitpre_gg_all[used_id_gg_plot], dr0_all_testing_gg_all[used_id_gg_plot], xlim = c(-10, 15), ylim = c(-10, 15), main=round(R2(exp_fitpre_gg_all[used_id_gg_plot], dr0_all_testing_gg_all[used_id_gg_plot]),3))
abline(0,1, col='red')
}
dev.off()
png('test.GG.split.P.LM.png', width=2000, height=500)
par(mfrow=c(1,4)) 
for (i in c(1:4)){
used_id_gg_plot = gene_exp_genegroup_vec_testing_new==i
plot(exp_fitpre_gg_all_P[used_id_gg_plot], dr0_all_testing_gg_all[used_id_gg_plot], xlim = c(-10, 15), ylim = c(-10, 15), main=round(R2(exp_fitpre_gg_all_P[used_id_gg_plot], dr0_all_testing_gg_all[used_id_gg_plot]),3))
abline(0,1, col='red')
}
dev.off()
png('test.GG.split.D.LM.png', width=2000, height=500)
par(mfrow=c(1,4)) 
for (i in c(1:4)){
used_id_gg_plot = gene_exp_genegroup_vec_testing_new==i
plot(exp_fitpre_gg_all_D[used_id_gg_plot], dr0_all_testing_gg_all[used_id_gg_plot], xlim = c(-10, 15), ylim = c(-10, 15), main=round(R2(exp_fitpre_gg_all_D[used_id_gg_plot], dr0_all_testing_gg_all[used_id_gg_plot]),3))
abline(0,1, col='red')
}
dev.off()




### DP
library(mclust)
fit_vs_obs = (cbind(exp_fitpre_gg_all_training, dr0_all_testing_gg_all_training))
BIC <- mclustBIC(fit_vs_obs, modelNames=c('EVV'))
summary(BIC)
mod1 <- Mclust(fit_vs_obs, x = BIC)
summary(mod1, parameters = TRUE)


png('test.GG.training.png', width=1000, height=500)
par(mfrow=c(1,2)) 
plot(exp_fitpre_gg_all_training, dr0_all_testing_gg_all_training, xlim = c(-10, 15), ylim = c(-10, 15), main=round(R2(exp_fitpre_gg_all_training, dr0_all_testing_gg_all_training),3))
abline(0,1, col='red')
plot(mod1, what = "classification", xlim = c(-10, 15), ylim = c(-10, 15))
dev.off()
png('test.GG.split.training.LM.png', width=2000, height=500)
par(mfrow=c(1,4)) 
for (i in c(1:4)){
used_id_gg_plot = gene_exp_genegroup_vec_new==i
plot(exp_fitpre_gg_all_training[used_id_gg_plot], dr0_all_testing_gg_all_training[used_id_gg_plot], xlim = c(-10, 15), ylim = c(-10, 15), main=round(R2(exp_fitpre_gg_all_training[used_id_gg_plot], dr0_all_testing_gg_all_training[used_id_gg_plot]),3))
abline(0,1, col='red')
}
dev.off()




### All GG
exp_fitpre_gg_all = c()
dr0_all_testing_gg_all = c()
exp_fitpre_gg_all_training2 = c()
dr0_all_testing_gg_all_training2 = c()
gene_exp_genegroup_vec1 = mod1$classification
adjustedRandIndex(gene_exp_genegroup_vec1, gene_exp_genegroup_vec_new)

new_cluster = apply(cbind(gene_exp_genegroup_vec1, gene_exp_genegroup_vec_new), 1, function(x) paste(x[1], x[2], collapse=''))
XXX = rownames(table(new_cluster))[table(new_cluster)<100]
toNum = cbind(rownames(table(new_cluster))[table(new_cluster)>=100])
toNum = cbind(toNum, seq(1, dim(toNum)[1]))
toNum = rbind(toNum, cbind(XXX, rep(dim(toNum)[1]+1, length(XXX)) ))

new_cluster_NUM = apply(cbind(new_cluster), 1, function(x) as.numeric(toNum[toNum[,1]==x,2]) )


for (gg in c(1:length(unique(new_cluster_NUM)))){
	print(gg)
	print(sum(new_cluster_NUM==gg))
dr0_all_gg = dr0_all_new[new_cluster_NUM==gg,]
fit_gg = lm(Gene_Group ~ .-1, data=dr0_all_gg[,])
exp_fitpre_gg_all_training2 = c(exp_fitpre_gg_all_training2, fit_gg$fitted.values)
dr0_all_testing_gg_all_training2 = c(dr0_all_testing_gg_all_training2, dr0_all_gg[,55])
fit_R2_gg = summary(fit_gg)$r.squared
print(fit_R2_gg)
dr0_all_testing_gg = dr0_all_testing[new_cluster_NUM[1:3414]==gg,]
print(dim(dr0_all_testing_gg))
exp_fitpre_gg = predict(fit_gg, dr0_all_testing_gg[,])
fit_R2pre_gg = R2(exp_fitpre_gg, dr0_all_testing_gg[,55])
exp_fitpre_gg_all = c(exp_fitpre_gg_all, exp_fitpre_gg)
dr0_all_testing_gg_all = c(dr0_all_testing_gg_all, dr0_all_testing_gg[,55])
print(fit_R2pre_gg)
}
fitpre_LM_R2_gg1 = R2(exp_fitpre_gg_all, dr0_all_testing_gg_all)
fitpre_LM_R2_gg1
0.7571352


fit_vs_obs2 = (cbind(exp_fitpre_gg_all_training2, dr0_all_testing_gg_all_training2))
BIC <- mclustBIC(fit_vs_obs2, modelNames=c('EVV'))
summary(BIC)
mod2 <- Mclust(fit_vs_obs2, x = BIC)
summary(mod2, parameters = TRUE)

png('test.GG.training.2.png', width=1000, height=500)
par(mfrow=c(1,2)) 
plot(exp_fitpre_gg_all_training2, dr0_all_testing_gg_all_training2, xlim = c(-10, 15), ylim = c(-10, 15), main=round(R2(exp_fitpre_gg_all_training2, dr0_all_testing_gg_all_training2),3))
abline(0,1, col='red')
plot(mod2, what = "classification", xlim = c(-10, 15), ylim = c(-10, 15))
dev.off()





### All GG _P
exp_fitpre_gg_all_P = c()
dr0_all_testing_gg_all_P = c()
for (gg in c(1:4)){
	print(gg)
	print(sum(gene_exp_genegroup_vec==gg))
dr0_all_gg = dr0_all[gene_exp_genegroup_vec==gg,]
fit_gg_P = lm(Gene_Group ~ .-1, data=dr0_all_gg[,28:55])
fit_R2_gg_P = summary(fit_gg_P)$r.squared
print(fit_R2_gg_P)
dr0_all_testing_gg_P = dr0_all_testing[gene_exp_genegroup_test==gg,]
exp_fitpre_gg_P = predict(fit_gg_P, dr0_all_testing_gg_P[,28:55])
fit_R2pre_gg_P = R2(exp_fitpre_gg_P, dr0_all_testing_gg_P[,55])
exp_fitpre_gg_all_P = c(exp_fitpre_gg_all_P, exp_fitpre_gg_P)
dr0_all_testing_gg_all_P = c(dr0_all_testing_gg_all_P, dr0_all_testing_gg_P[,55])
print(fit_R2pre_gg_P)
}
fitpre_LM_R2_gg_P = R2(exp_fitpre_gg_all_P, dr0_all_testing_gg_all_P)
fitpre_LM_R2_gg_P
0.7169234
### All GG _D
exp_fitpre_gg_all_D = c()
dr0_all_testing_gg_all_D = c()
for (gg in c(1:4)){
	print(gg)
	print(sum(gene_exp_genegroup_vec==gg))
dr0_all_gg = dr0_all[gene_exp_genegroup_vec==gg,]
fit_gg_D = lm(Gene_Group ~ .-1, data=dr0_all_gg[,c(1:27,55)])
fit_R2_gg_D = summary(fit_gg_D)$r.squared
print(fit_R2_gg_D)
dr0_all_testing_gg_D = dr0_all_testing[gene_exp_genegroup_test==gg,]
exp_fitpre_gg_D = predict(fit_gg_D, dr0_all_testing_gg_D[,c(1:27,55)])
fit_R2pre_gg_D = R2(exp_fitpre_gg_D, dr0_all_testing_gg_D[,55])
exp_fitpre_gg_all_D = c(exp_fitpre_gg_all_D, exp_fitpre_gg_D)
dr0_all_testing_gg_all_D = c(dr0_all_testing_gg_all_D, dr0_all_testing_gg_D[,55])
print(fit_R2pre_gg_D)
}
fitpre_LM_R2_gg_D = R2(exp_fitpre_gg_all_D, dr0_all_testing_gg_all_D)
fitpre_LM_R2_gg_D
0.7450707













for (i in c(1:4)){
	print(i)
	dr0 = c()
	gene_exp_tmp = c()
	for (j in c(1:11)){
		print(j)
		gene_exp_tmp_tmp = c()
		for (k in chr_num){
		print(k)
		load(paste('../vision_rna_tss2k_ccreunit.chr', k, '.', i, '.', j, '.Rdat', sep=''))
		used_id_ct = c(1:(length(rt$y)/12)) + (j-1)*(length(rt$y)/12)
		dr0 = rbind(dr0, cbind(rt$x[used_id_ct,], rt$x0[used_id_ct,], rt$y[used_id_ct]))
		label_all_genegroup0 = c(label_all_genegroup0, rep(i, length(rt$y[used_id_ct])))
		gene_exp_tmp_tmp = c(gene_exp_tmp_tmp, rt$y[used_id_ct])
		}
		gene_exp_tmp = cbind(gene_exp_tmp, gene_exp_tmp_tmp)
	}
	gene_exp = rbind(gene_exp, gene_exp_tmp)
	fit_lm = lm(dr0[,55]~dr0[,-55]-1)
	eRP = fit_lm$coefficients
	eRP_mat_D = cbind(eRP_mat_D, eRP[1:27])
	eRP_mat_P = cbind(eRP_mat_P, eRP[28:54])
	dr_kmeans_plot_cor = cor(dr0, method='pearson')
	dr_kmeans_plot_cor_all = cbind(dr_kmeans_plot_cor_all, dr_kmeans_plot_cor[-55,55])
	dr0_all = rbind(dr0_all, dr0)
}
eRP_mat = cbind(eRP_mat_D, eRP_mat_P)

m = apply(gene_exp, 1, mean)
s = apply(gene_exp, 1, sd)

tt3 = which(m>-4 & s>2)
tt4 = which(m>-4 & s<=2)
tt2 = which(m<=-4 & s>2)
tt1 = which(m<=-4 & s<=2)

label_all_genegroup1 = rep(0, length(m))
label_all_genegroup1[tt3] = 3
label_all_genegroup1[tt4] = 4
label_all_genegroup1[tt2] = 2
label_all_genegroup1[tt1] = 1

label_all_genegroup = rep(label_all_genegroup1, 11)

set.seed(2018)
kmeans_fit = kmeans(gene_exp, centers=20)
breaksList_scale = seq(0, max(gene_exp), by = 0.01)
print(breaksList_scale)
my_colorbar_scale=colorRampPalette(c('white', 'red'))(n = length(breaksList_scale))


png('kmeans.rna.tpm.png')
pheatmap(gene_exp[order(kmeans_fit$cluster),], color=my_colorbar_scale, cluster_rows = F, cluster_cols = F,annotation_names_row=FALSE,annotation_names_col=FALSE,show_rownames=F,show_colnames=F)
dev.off()


dr0_all_testing = c()
#chr_num = c(1:19, 'X')
chr_num_testing = c(1)
for (i in c(1:4)){
	print(i)
	dr0 = c()
	for (j in c(12)){
		print(j)
		for (k in chr_num_testing){
		print(k)
		load(paste('../vision_rna_tss2k_ccreunit.chr', k, '.', i, '.', j, '.Rdat', sep=''))
		used_id_ct = c(1:(length(rt$y)/12)) + (j-1)*(length(rt$y)/12)
		dr0 = rbind(dr0, cbind(rt$x[used_id_ct,], rt$x0[used_id_ct,], rt$y[used_id_ct]))
		}
	}
	dr0_all_testing = rbind(dr0_all_testing, dr0)
}




colnames_state = c()
for (i in 0:26){
	colnames_state = c(colnames_state, paste('D',i,sep=''))
}
for (i in 0:26){
	colnames_state = c(colnames_state, paste('P',i,sep=''))
}
colnames_state = c(colnames_state, 'Gene_Group')


colnames(dr0_all) = colnames_state
colnames(dr0_all_testing) = colnames_state
dr0_all = as.data.frame(dr0_all)
dr0_all_testing = as.data.frame(dr0_all_testing)

set.seed(2018)
used_id_chr1_2 =  sample(dim(dr0_all)[1], 40000)
used_id_chr3_4 =  sample(dim(dr0_all_testing)[1], 40000)
used_id_chr1_2_nn =  sample(dim(dr0_all)[1], 1000)

### All
dr0_all_training_all = c()
dr0_all_testing_all = c()
for (gg in c(1:4)){
dr0_all_gg = dr0_all[label_all_genegroup==gg,]
dr0_all_training_all = rbind(dr0_all_training_all, data=dr0_all_gg[,])
##
dr0_all_testing_gg = dr0_all_testing[label_all_genegroup1==gg,]
dr0_all_testing_all = rbind(dr0_all_testing_all, dr0_all_testing_gg[,])
}

fit = lm(Gene_Group ~ .-1, data=dr0_all_training_all)
fit_R2 = summary(fit)$r.squared
fit_R2
exp_fitpre = predict(fit, dr0_all_testing_all[,-55])
fit_R2pre = R2(exp_fitpre, dr0_all_testing_all[,55])
fit_R2pre
0.4876924
0.5108249
0.5171285
### Proximal
fit_P = lm(Gene_Group ~ .-1, data=dr0_all_training_all[,28:55])
fit_P_R2 = summary(fit_P)$r.squared
fit_P_R2
exp_fitpre_P = predict(fit_P, dr0_all_testing_all[,28:54])
fit_R2pre_P = R2(exp_fitpre_P, dr0_all_testing_all[,55])
fit_R2pre_P
0.334773
0.4138006
0.4726581
### Distal
fit_D = lm(Gene_Group ~ .-1, data=dr0_all_training_all[,c(1:27,55)])
fit_D_R2 = summary(fit_D)$r.squared
fit_D_R2
exp_fitpre_D = predict(fit_D, dr0_all_testing_all[,c(1:27,55)])
fit_R2pre_D = R2(exp_fitpre_D, dr0_all_testing_all[,55])
fit_R2pre_D
0.3228512
0.3278979
0.2542589

png('test.LM.png', width=1500, height=500)
par(mfrow=c(1,3)) 
plot(exp_fitpre, dr0_all_testing_all[,55], xlim = c(-10, 15), ylim = c(-10, 15), main=round(R2(exp_fitpre, dr0_all_testing_all[,55]), 3))
abline(0,1, col='red')
plot(exp_fitpre_P, dr0_all_testing_all[,55], xlim = c(-10, 15), ylim = c(-10, 15), main=round(R2(exp_fitpre_P, dr0_all_testing_all[,55]), 3))
abline(0,1, col='red')
plot(exp_fitpre_D, dr0_all_testing_all[,55], xlim = c(-10, 15), ylim = c(-10, 15), main=round(R2(exp_fitpre_D, dr0_all_testing_all[,55]), 3))
abline(0,1, col='red')
dev.off()


### All GG
exp_fitpre_gg_all = c()
dr0_all_testing_gg_all = c()
exp_fitpre_gg_all_training = c()
dr0_all_testing_gg_all_training = c()
for (gg in c(1:4)){
	print(gg)
	print(sum(label_all_genegroup==gg))
dr0_all_gg = dr0_all[label_all_genegroup==gg,]
fit_gg = lm(Gene_Group ~ .-1, data=dr0_all_gg[,])
exp_fitpre_gg_all_training = rbind(exp_fitpre_gg_all_training, fit_gg$fitted.values)
dr0_all_testing_gg_all_training = rbind(dr0_all_testing_gg_all_training, dr0_all_gg[,55])
fit_R2_gg = summary(fit_gg)$r.squared
print(fit_R2_gg)
dr0_all_testing_gg = dr0_all_testing[label_all_genegroup1==gg,]
exp_fitpre_gg = predict(fit_gg, dr0_all_testing_gg[,])
fit_R2pre_gg = R2(exp_fitpre_gg, dr0_all_testing_gg[,55])
exp_fitpre_gg_all = c(exp_fitpre_gg_all, exp_fitpre_gg)
dr0_all_testing_gg_all = c(dr0_all_testing_gg_all, dr0_all_testing_gg[,55])
print(fit_R2pre_gg)
}
fitpre_LM_R2_gg = R2(exp_fitpre_gg_all, dr0_all_testing_gg_all)
fitpre_LM_R2_gg
0.7496579
0.7519384
0.8226268
0.8233064

### All Proximal
exp_fitpre_gg_all_P = c()
dr0_all_testing_gg_all_P = c()
for (gg in c(1:4)){
#gg = 1
dr0_all_gg = dr0_all[label_all_genegroup==gg,]
fit_gg = lm(Gene_Group ~ .-1, data=dr0_all_gg[,c(28:54,55)])
fit_R2_gg = summary(fit_gg)$r.squared
print(fit_R2_gg)
dr0_all_testing_gg = dr0_all_testing[label_all_genegroup1==gg,]
exp_fitpre_gg = predict(fit_gg, dr0_all_testing_gg[,c(28:54,55)])
fit_R2pre_gg_P = R2(exp_fitpre_gg, dr0_all_testing_gg[,55])
exp_fitpre_gg_all_P = c(exp_fitpre_gg_all_P, exp_fitpre_gg)
dr0_all_testing_gg_all_P = c(dr0_all_testing_gg_all_P, dr0_all_testing_gg[,55])
print(fit_R2pre_gg_P)
}
fitpre_LM_R2_gg_P = R2(exp_fitpre_gg_all_P, dr0_all_testing_gg_all_P)
fitpre_LM_R2_gg_P
0.6351115
0.641122
0.7259135
0.7273904

### All Distal
exp_fitpre_gg_all_D = c()
dr0_all_testing_gg_all_D = c()
for (gg in c(1:4)){
#gg = 1
dr0_all_gg = dr0_all[label_all_genegroup==gg,]
set.seed(2018)
fit_gg = lm(Gene_Group ~ .-1, data=dr0_all_gg[,c(1:27,55)])
fit_R2_gg = summary(fit_gg)$r.squared
print(fit_R2_gg)
dr0_all_testing_gg = dr0_all_testing[label_all_genegroup1==gg,]
exp_fitpre_gg = predict(fit_gg, dr0_all_testing_gg[,c(1:27,55)])
fit_R2pre_gg_D = R2(exp_fitpre_gg, dr0_all_testing_gg[,55])
exp_fitpre_gg_all_D = c(exp_fitpre_gg_all_D, exp_fitpre_gg)
dr0_all_testing_gg_all_D = c(dr0_all_testing_gg_all_D, dr0_all_testing_gg[,55])
print(fit_R2pre_gg_D)
}
fitpre_LM_R2_gg_D = R2(exp_fitpre_gg_all_D, dr0_all_testing_gg_all_D)
fitpre_LM_R2_gg_D
0.7177247
0.7182962
0.8162241
0.8163163


png('test.GG.training.png', width=1500, height=500)
par(mfrow=c(1,3)) 
plot(exp_fitpre_gg_all_training, dr0_all_testing_gg_all_training, xlim = c(-10, 15), ylim = c(-10, 15), main=round(R2(exp_fitpre_gg_all, dr0_all_testing_gg_all),3))
abline(0,1, col='red')
plot(exp_fitpre_gg_all_P, dr0_all_testing_gg_all_P, xlim = c(-10, 15), ylim = c(-10, 15), main=round(R2(exp_fitpre_gg_all_P, dr0_all_testing_gg_all_P),3))
abline(0,1, col='red')
plot(exp_fitpre_gg_all_D, dr0_all_testing_gg_all_D, xlim = c(-10, 15), ylim = c(-10, 15), main=round(R2(exp_fitpre_gg_all_D, dr0_all_testing_gg_all_D),3))
abline(0,1, col='red')
dev.off()

png('test.GG.png', width=1500, height=500)
par(mfrow=c(1,3)) 
plot(exp_fitpre_gg_all, dr0_all_testing_gg_all, xlim = c(-10, 15), ylim = c(-10, 15), main=round(R2(exp_fitpre_gg_all, dr0_all_testing_gg_all),3))
abline(0,1, col='red')
plot(exp_fitpre_gg_all_P, dr0_all_testing_gg_all_P, xlim = c(-10, 15), ylim = c(-10, 15), main=round(R2(exp_fitpre_gg_all_P, dr0_all_testing_gg_all_P),3))
abline(0,1, col='red')
plot(exp_fitpre_gg_all_D, dr0_all_testing_gg_all_D, xlim = c(-10, 15), ylim = c(-10, 15), main=round(R2(exp_fitpre_gg_all_D, dr0_all_testing_gg_all_D),3))
abline(0,1, col='red')
dev.off()


### All GG
exp_fitpre_gg_all_lm = c()
dr0_all_testing_gg_all_lm = c()
for (gg in c(1:4)){
#gg = 1
dr0_all_gg = dr0_all[label_all_genegroup==gg,]
set.seed(2018)
used_id_gg = sample(dim(dr0_all_gg)[1], 10000)
fit_gg = lm(Gene_Group ~ .-1, data=dr0_all_gg[used_id_gg,])
fit_R2_gg = summary(fit_gg)$r.squared
print(fit_R2_gg)
dr0_all_testing_gg = dr0_all_gg[-used_id_gg,]
used_id_chr3_4_gg = sample(dim(dr0_all_testing_gg)[1], 10000)
exp_fitpre_gg = predict(fit_gg, dr0_all_testing_gg[used_id_chr3_4_gg,])
fit_R2pre_gg_lm = R2(exp_fitpre_gg, dr0_all_testing_gg[used_id_chr3_4_gg,55])
exp_fitpre_gg_all_lm = c(exp_fitpre_gg_all_lm, exp_fitpre_gg)
dr0_all_testing_gg_all_lm = c(dr0_all_testing_gg_all_lm, dr0_all_testing_gg[used_id_chr3_4_gg,55])
print(fit_R2pre_gg_lm)
}
fitpre_LM_R2_gg_lm = R2(exp_fitpre_gg_all_lm, dr0_all_testing_gg_all_lm)
fitpre_LM_R2_gg_lm
0.7496579

### All Proximal
exp_fitpre_gg_all_P_lm = c()
dr0_all_testing_gg_all_P_lm = c()
for (gg in c(1:4)){
#gg = 1
dr0_all_gg = dr0_all[label_all_genegroup==gg,]
set.seed(2018)
used_id_gg = sample(dim(dr0_all_gg)[1], 10000)
fit_gg = lm(Gene_Group ~ .-1, data=dr0_all_gg[used_id_gg,c(28:54,55)])
fit_R2_gg = summary(fit_gg)$r.squared
print(fit_R2_gg)
dr0_all_testing_gg = dr0_all_gg[-used_id_gg,]
used_id_chr3_4_gg = sample(dim(dr0_all_testing_gg)[1], 10000)
exp_fitpre_gg = predict(fit_gg, dr0_all_testing_gg[used_id_chr3_4_gg,c(28:54,55)])
fit_R2pre_gg_P_lm = R2(exp_fitpre_gg, dr0_all_testing_gg[used_id_chr3_4_gg,55])
exp_fitpre_gg_all_P_lm = c(exp_fitpre_gg_all_P_lm, exp_fitpre_gg)
dr0_all_testing_gg_all_P_lm = c(dr0_all_testing_gg_all_P_lm, dr0_all_testing_gg[used_id_chr3_4_gg,55])
print(fit_R2pre_gg_P_lm)
}
fitpre_LM_R2_gg_P_lm = R2(exp_fitpre_gg_all_P_lm, dr0_all_testing_gg_all_P_lm)
fitpre_LM_R2_gg_P_lm
0.6351115

### All Distal
exp_fitpre_gg_all_D_lm = c()
dr0_all_testing_gg_all_D_lm = c()
for (gg in c(1:4)){
#gg = 1
dr0_all_gg = dr0_all[label_all_genegroup==gg,]
set.seed(2018)
used_id_gg = sample(dim(dr0_all_gg)[1], 10000)
fit_gg = lm(Gene_Group ~ .-1, data=dr0_all_gg[used_id_gg,c(1:27,55)])
fit_R2_gg = summary(fit_gg)$r.squared
print(fit_R2_gg)
dr0_all_testing_gg = dr0_all_gg[-used_id_gg,]
used_id_chr3_4_gg = sample(dim(dr0_all_testing_gg)[1], 10000)
exp_fitpre_gg = predict(fit_gg, dr0_all_testing_gg[used_id_chr3_4_gg,c(1:27,55)])
fit_R2pre_gg_D_lm = R2(exp_fitpre_gg, dr0_all_testing_gg[used_id_chr3_4_gg,55])
exp_fitpre_gg_all_D_lm = c(exp_fitpre_gg_all_D_lm, exp_fitpre_gg)
dr0_all_testing_gg_all_D_lm = c(dr0_all_testing_gg_all_D_lm, dr0_all_testing_gg[used_id_chr3_4_gg,55])
print(fit_R2pre_gg_D_lm)
}
fitpre_LM_R2_gg_D_lm = R2(exp_fitpre_gg_all_D_lm, dr0_all_testing_gg_all_D_lm)
fitpre_LM_R2_gg_D_lm
0.7177247


png('test.GG.split.LM.png', width=2000, height=500)
par(mfrow=c(1,4)) 
used_id_gg_plot = c(1:10000)
plot(exp_fitpre_gg_all_lm[used_id_gg_plot], dr0_all_testing_gg_all_lm[used_id_gg_plot], xlim = c(-10, 15), ylim = c(-10, 15), main=round(R2(exp_fitpre_gg_all_lm[used_id_gg_plot], dr0_all_testing_gg_all_lm[used_id_gg_plot]),3))
abline(0,1, col='red')
used_id_gg_plot = c(10001:20000)
plot(exp_fitpre_gg_all_lm[used_id_gg_plot], dr0_all_testing_gg_all_lm[used_id_gg_plot], xlim = c(-10, 15), ylim = c(-10, 15), main=round(R2(exp_fitpre_gg_all_lm[used_id_gg_plot], dr0_all_testing_gg_all_lm[used_id_gg_plot]),3))
abline(0,1, col='red')
used_id_gg_plot = c(20001:30000)
plot(exp_fitpre_gg_all_lm[used_id_gg_plot], dr0_all_testing_gg_all_lm[used_id_gg_plot], xlim = c(-10, 15), ylim = c(-10, 15), main=round(R2(exp_fitpre_gg_all_lm[used_id_gg_plot], dr0_all_testing_gg_all_lm[used_id_gg_plot]),3))
abline(0,1, col='red')
used_id_gg_plot = c(30001:40000)
plot(exp_fitpre_gg_all_lm[used_id_gg_plot], dr0_all_testing_gg_all_lm[used_id_gg_plot], xlim = c(-10, 15), ylim = c(-10, 15), main=round(R2(exp_fitpre_gg_all_lm[used_id_gg_plot], dr0_all_testing_gg_all_lm[used_id_gg_plot]),3))
abline(0,1, col='red')
dev.off()


png('test.GG.split.D.LM.png', width=2000, height=500)
par(mfrow=c(1,4)) 
used_id_gg_plot = c(1:10000)
plot(exp_fitpre_gg_all_D_lm[used_id_gg_plot], dr0_all_testing_gg_all_D_lm[used_id_gg_plot], xlim = c(-10, 15), ylim = c(-10, 15), main=round(R2(exp_fitpre_gg_all_D_lm[used_id_gg_plot], dr0_all_testing_gg_all_D_lm[used_id_gg_plot]),3))
abline(0,1, col='red')
used_id_gg_plot = c(10001:20000)
plot(exp_fitpre_gg_all_D_lm[used_id_gg_plot], dr0_all_testing_gg_all_D_lm[used_id_gg_plot], xlim = c(-10, 15), ylim = c(-10, 15), main=round(R2(exp_fitpre_gg_all_D_lm[used_id_gg_plot], dr0_all_testing_gg_all_D_lm[used_id_gg_plot]),3))
abline(0,1, col='red')
used_id_gg_plot = c(20001:30000)
plot(exp_fitpre_gg_all_D_lm[used_id_gg_plot], dr0_all_testing_gg_all_D_lm[used_id_gg_plot], xlim = c(-10, 15), ylim = c(-10, 15), main=round(R2(exp_fitpre_gg_all_D_lm[used_id_gg_plot], dr0_all_testing_gg_all_D_lm[used_id_gg_plot]),3))
abline(0,1, col='red')
used_id_gg_plot = c(30001:40000)
plot(exp_fitpre_gg_all_D_lm[used_id_gg_plot], dr0_all_testing_gg_all_D_lm[used_id_gg_plot], xlim = c(-10, 15), ylim = c(-10, 15), main=round(R2(exp_fitpre_gg_all_D_lm[used_id_gg_plot], dr0_all_testing_gg_all_D_lm[used_id_gg_plot]),3))
abline(0,1, col='red')
dev.off()


png('test.GG.split.P.LM.png', width=2000, height=500)
par(mfrow=c(1,4)) 
used_id_gg_plot = c(1:10000)
plot(exp_fitpre_gg_all_P_lm[used_id_gg_plot], dr0_all_testing_gg_all_P_lm[used_id_gg_plot], xlim = c(-10, 15), ylim = c(-10, 15), main=round(R2(exp_fitpre_gg_all_P_lm[used_id_gg_plot], dr0_all_testing_gg_all_P_lm[used_id_gg_plot]),3))
abline(0,1, col='red')
used_id_gg_plot = c(10001:20000)
plot(exp_fitpre_gg_all_P_lm[used_id_gg_plot], dr0_all_testing_gg_all_P_lm[used_id_gg_plot], xlim = c(-10, 15), ylim = c(-10, 15), main=round(R2(exp_fitpre_gg_all_P_lm[used_id_gg_plot], dr0_all_testing_gg_all_P_lm[used_id_gg_plot]),3))
abline(0,1, col='red')
used_id_gg_plot = c(20001:30000)
plot(exp_fitpre_gg_all_P_lm[used_id_gg_plot], dr0_all_testing_gg_all_P_lm[used_id_gg_plot], xlim = c(-10, 15), ylim = c(-10, 15), main=round(R2(exp_fitpre_gg_all_P_lm[used_id_gg_plot], dr0_all_testing_gg_all_P_lm[used_id_gg_plot]),3))
abline(0,1, col='red')
used_id_gg_plot = c(30001:40000)
plot(exp_fitpre_gg_all_P_lm[used_id_gg_plot], dr0_all_testing_gg_all_P_lm[used_id_gg_plot], xlim = c(-10, 15), ylim = c(-10, 15), main=round(R2(exp_fitpre_gg_all_P_lm[used_id_gg_plot], dr0_all_testing_gg_all_P_lm[used_id_gg_plot]),3))
abline(0,1, col='red')
dev.off()

png('test.GG.LM.png', width=1500, height=500)
par(mfrow=c(1,3)) 
plot(exp_fitpre_gg_all_lm, dr0_all_testing_gg_all_lm, xlim = c(-10, 15), ylim = c(-10, 15), main=round(R2(exp_fitpre_gg_all_lm, dr0_all_testing_gg_all_lm),3))
abline(0,1, col='red')
plot(exp_fitpre_gg_all_P_lm, dr0_all_testing_gg_all_P_lm, xlim = c(-10, 15), ylim = c(-10, 15), main=round(R2(exp_fitpre_gg_all_P_lm, dr0_all_testing_gg_all_P_lm),3))
abline(0,1, col='red')
plot(exp_fitpre_gg_all_D_lm, dr0_all_testing_gg_all_D_lm, xlim = c(-10, 15), ylim = c(-10, 15), main=round(R2(exp_fitpre_gg_all_D_lm, dr0_all_testing_gg_all_D_lm),3))
abline(0,1, col='red')
dev.off()


