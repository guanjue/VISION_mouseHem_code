library(randomForest)
library(pheatmap)
library(e1071)
library(neuralnet)
library(LSD)

###########
R2 = function(d_IP_ts, qPCR_expand){
        r2 = 1-sum((qPCR_expand-d_IP_ts)^2) / sum((qPCR_expand-mean(qPCR_expand))^2)
        return(r2)
}
###########

sigmat = read.table('~/group/projects/vision/pknorm_2_16lim_ref1mo_sample1mo_0424/pknorm_2_16lim_ref1mo_0424_lesshet.para0')
sigmat1 = as.matrix(sigmat[,2:9]/sigmat[,1])
print(head(sigmat1))
fit=hclust(dist(sigmat1),method="ward.D2")


dr_kmeans_plot_cor_all = c()
dr0_all = c()
gene_exp = c()
label_all_genegroup0 = c()
eRP_mat_D = c()
eRP_mat_P = c()
#chr_num = c(1:19, 'X')
chr_num = c(1)
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
		dr0 = rbind(dr0, cbind(rt$x, rt$x0, rt$y))
		label_all_genegroup0 = c(label_all_genegroup0, rep(i, length(rt$y)))
		gene_exp_tmp_tmp = c(gene_exp_tmp_tmp, rt$y)
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

label_all_genegroup = label_all_genegroup0
label_all_genegroup[tt3] = 3
label_all_genegroup[tt4] = 4
label_all_genegroup[tt2] = 2
label_all_genegroup[tt1] = 1


set.seed(2018)
kmeans_fit = kmeans(gene_exp, centers=20)
my_colorbar_scale=colorRampPalette(c('white', 'red'))(n = length(breaksList_scale))

png('kmeans.rna.tpm.png')
pheatmap(gene_exp[order(kmeans_fit$cluster),], color=my_colorbar_scale, cluster_rows = F, cluster_cols = F,annotation_names_row=FALSE,annotation_names_col=FALSE,show_rownames=F,show_colnames=F)
dev.off()

dr0_all_testing = c()
#chr_num = c(1:19, 'X')
chr_num_testing = c(2)
for (i in c(1:4)){
	print(i)
	dr0 = c()
	for (j in c(12)){
		print(j)
		for (k in chr_num_testing){
		print(k)
		load(paste('../vision_rna_tss2k_ccreunit.chr', k, '.', i, '.', j, '.Rdat', sep=''))
		dr0 = rbind(dr0, cbind(rt$x, rt$x0, rt$y))
		}
	}
	dr0_all_testing = rbind(dr0_all_testing, dr0)
}



group_cluster = c()
chr_num_all = c()
for (i in c(1:4)){
	print(i)
	dr0 = c()
	for (j in c(1:12)){
		print(j)
		chr_num_tmp = c()
		for (k in chr_num){
		#print(k)
		load(paste('../vision_rna_tss2k_ccreunit.chr', k, '.', i, '.', j, '.Rdat', sep=''))
		tmp = length(rt$y)
		#print(tmp)
		chr_num_tmp = c(chr_num_tmp, tmp)
		}
		#print(chr_num)
		chr_num_all = rbind(chr_num_all, chr_num_tmp)
	}
}


group_num_all = colSums(matrix(apply(as.matrix(chr_num_all), 1, function(x) sum(as.numeric(x))), nrow=12))
label = c()
for ( i in c(1:4)){
	label = c(label, rep(i, group_num_all[i]))
}
cater_sig_label = cbind(dr0_all[,1:54], label)
cater = cor(cater_sig_label)
cater_mat = cbind(cater[1:27,55], cater[28:54,55])


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
set.seed(2018)
used_id_gg = sample(dim(dr0_all_gg)[1], 10000)
dr0_all_training_all = rbind(dr0_all_training_all, data=dr0_all_gg[used_id_gg,])
dr0_all_testing_gg = dr0_all_gg[-used_id_gg,]
used_id_chr3_4_gg = sample(dim(dr0_all_testing_gg)[1], 10000)
dr0_all_testing_all = rbind(dr0_all_testing_all, dr0_all_testing_gg[used_id_chr3_4_gg,])
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
for (gg in c(1:4)){
	print(gg)
	print(sum(label_all_genegroup==gg))
dr0_all_gg = dr0_all[label_all_genegroup==gg,]
set.seed(2018)
used_id_gg = sample(dim(dr0_all_gg)[1], 10000)
fit_gg = lm(Gene_Group ~ .-1, data=dr0_all_gg[used_id_gg,])
fit_R2_gg = summary(fit_gg)$r.squared
print(fit_R2_gg)
dr0_all_testing_gg = dr0_all_gg[-used_id_gg,]
used_id_chr3_4_gg = sample(dim(dr0_all_testing_gg)[1], 10000)
exp_fitpre_gg = predict(fit_gg, dr0_all_testing_gg[used_id_chr3_4_gg,])
fit_R2pre_gg = R2(exp_fitpre_gg, dr0_all_testing_gg[used_id_chr3_4_gg,55])
exp_fitpre_gg_all = c(exp_fitpre_gg_all, exp_fitpre_gg)
dr0_all_testing_gg_all = c(dr0_all_testing_gg_all, dr0_all_testing_gg[used_id_chr3_4_gg,55])
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
set.seed(2018)
used_id_gg = sample(dim(dr0_all_gg)[1], 10000)
fit_gg = lm(Gene_Group ~ .-1, data=dr0_all_gg[used_id_gg,c(28:54,55)])
fit_R2_gg = summary(fit_gg)$r.squared
print(fit_R2_gg)
dr0_all_testing_gg = dr0_all_gg[-used_id_gg,]
used_id_chr3_4_gg = sample(dim(dr0_all_testing_gg)[1], 10000)
exp_fitpre_gg = predict(fit_gg, dr0_all_testing_gg[used_id_chr3_4_gg,c(28:54,55)])
fit_R2pre_gg_P = R2(exp_fitpre_gg, dr0_all_testing_gg[used_id_chr3_4_gg,55])
exp_fitpre_gg_all_P = c(exp_fitpre_gg_all_P, exp_fitpre_gg)
dr0_all_testing_gg_all_P = c(dr0_all_testing_gg_all_P, dr0_all_testing_gg[used_id_chr3_4_gg,55])
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
used_id_gg = sample(dim(dr0_all_gg)[1], 10000)
fit_gg = lm(Gene_Group ~ .-1, data=dr0_all_gg[used_id_gg,c(1:27,55)])
fit_R2_gg = summary(fit_gg)$r.squared
print(fit_R2_gg)
dr0_all_testing_gg = dr0_all_gg[-used_id_gg,]
used_id_chr3_4_gg = sample(dim(dr0_all_testing_gg)[1], 10000)
exp_fitpre_gg = predict(fit_gg, dr0_all_testing_gg[used_id_chr3_4_gg,c(1:27,55)])
fit_R2pre_gg_D = R2(exp_fitpre_gg, dr0_all_testing_gg[used_id_chr3_4_gg,55])
exp_fitpre_gg_all_D = c(exp_fitpre_gg_all_D, exp_fitpre_gg)
dr0_all_testing_gg_all_D = c(dr0_all_testing_gg_all_D, dr0_all_testing_gg[used_id_chr3_4_gg,55])
print(fit_R2pre_gg_D)
}
fitpre_LM_R2_gg_D = R2(exp_fitpre_gg_all_D, dr0_all_testing_gg_all_D)
fitpre_LM_R2_gg_D
0.7177247
0.7182962
0.8162241
0.8163163


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


