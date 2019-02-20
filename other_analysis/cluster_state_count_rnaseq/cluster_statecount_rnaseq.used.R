library(randomForest)
library(pheatmap)
library(e1071)
library(neuralnet)
library(LSD)
library(lcmm)

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
chr_num = c(1:19, 'X')
chr_num = c(1,3:19, 'X')
#chr_num = c(1)
eRP_mat_D = c()
eRP_mat_P = c()

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
fit = lm(Gene_Group ~ .-1, data=dr0_all[used_id_chr1_2,])
fit_R2 = summary(fit)$r.squared
fit_R2
exp_fitpre = predict(fit, as.data.frame(dr0_all_testing[used_id_chr3_4,-55]))
fit_R2pre = R2(exp_fitpre, dr0_all_testing[used_id_chr3_4,55])
fit_R2pre
0.5221678
0.5108249
0.5171285
### Proximal
fit_P = lm(Gene_Group ~ .-1, data=dr0_all[used_id_chr1_2,28:55])
fit_P_R2 = summary(fit_P)$r.squared
fit_P_R2
exp_fitpre_P = predict(fit_P, as.data.frame(dr0_all_testing[used_id_chr3_4,28:54]))
fit_R2pre_P = R2(exp_fitpre_P, dr0_all_testing[used_id_chr3_4,55])
fit_R2pre_P
0.4194863
0.4138006
0.4726581
### Distal
fit_D = lm(Gene_Group ~ .-1, data=dr0_all[used_id_chr1_2,c(1:27,55)])
fit_D_R2 = summary(fit_D)$r.squared
fit_D_R2
exp_fitpre_D = predict(fit_D, as.data.frame(dr0_all_testing[used_id_chr3_4,c(1:27,55)]))
fit_R2pre_D = R2(exp_fitpre_D, dr0_all_testing[used_id_chr3_4,55])
fit_R2pre_D
0.3352298
0.3278979
0.2542589

png('test.LM.png', width=1500, height=500)
par(mfrow=c(1,3)) 
plot(exp_fitpre, dr0_all_testing[used_id_chr3_4,55], xlim = c(-10, 15), ylim = c(-10, 15), main=round(R2(exp_fitpre, dr0_all_testing[used_id_chr3_4,55]), 3))
abline(0,1, col='red')
plot(exp_fitpre_P, dr0_all_testing[used_id_chr3_4,55], xlim = c(-10, 15), ylim = c(-10, 15), main=round(R2(exp_fitpre_P, dr0_all_testing[used_id_chr3_4,55]), 3))
abline(0,1, col='red')
plot(exp_fitpre_D, dr0_all_testing[used_id_chr3_4,55], xlim = c(-10, 15), ylim = c(-10, 15), main=round(R2(exp_fitpre_D, dr0_all_testing[used_id_chr3_4,55]), 3))
abline(0,1, col='red')
dev.off()

### All SVR
modelsvm = svm(Gene_Group ~ ., data=dr0_all[used_id_chr1_2,])
fitpre_svr = predict(modelsvm, as.data.frame(dr0_all_testing[used_id_chr3_4,-55]))
fitpre_svrR2 = R2(fitpre_svr, dr0_all_testing[used_id_chr3_4,55])
fitpre_svrR2
0.5514396
0.6016268
### Promixal SVR
modelsvm_P = svm(Gene_Group ~ ., data=dr0_all[used_id_chr1_2,28:55])
fitpre_P_svr = predict(modelsvm_P, as.data.frame(dr0_all_testing[used_id_chr3_4,28:55]))
fitpre_P_svrR2 = R2(fitpre_P_svr, dr0_all_testing[used_id_chr3_4,55])
fitpre_P_svrR2
0.4402826
0.4352918
### Distal SVR
modelsvm_D = svm(Gene_Group ~ ., data=dr0_all[used_id_chr1_2,c(1:27,55)])
fitpre_D_svr = predict(modelsvm_D, as.data.frame(dr0_all_testing[used_id_chr3_4,c(1:27,55)]))
fitpre_D_svrR2 = R2(fitpre_D_svr, dr0_all_testing[used_id_chr3_4,55])
fitpre_D_svrR2
0.450381
0.4535795

png('test.SVR.png', width=1500, height=500)
par(mfrow=c(1,3)) 
plot(fitpre_svr, dr0_all_testing[used_id_chr3_4,55], xlim = c(-10, 15), ylim = c(-10, 15), main=round(R2(fitpre_svr, dr0_all_testing[used_id_chr3_4,55]), 3))
abline(0,1, col='red')
plot(fitpre_P_svr, dr0_all_testing[used_id_chr3_4,55], xlim = c(-10, 15), ylim = c(-10, 15), main=round(R2(fitpre_P_svr, dr0_all_testing[used_id_chr3_4,55]), 3))
abline(0,1, col='red')
plot(fitpre_D_svr, dr0_all_testing[used_id_chr3_4,55], xlim = c(-10, 15), ylim = c(-10, 15), main=round(R2(fitpre_D_svr, dr0_all_testing[used_id_chr3_4,55]), 3))
abline(0,1, col='red')
dev.off()


### All GG
exp_fitpre_gg_all = c()
dr0_all_testing_gg_all = c()
for (gg in c(1:4)){
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


### All GG SVM
exp_fitpre_gg_all_svm = c()
dr0_all_testing_gg_all_svm = c()
for (gg in c(1:4)){
#gg = 1
dr0_all_gg = dr0_all[label_all_genegroup==gg,]
set.seed(2018)
used_id_gg = sample(dim(dr0_all_gg)[1], 10000)
fit_gg = svm(Gene_Group ~ .-1, data=dr0_all_gg[used_id_gg,])
fit_R2_gg = summary(fit_gg)$r.squared
print(fit_R2_gg)
dr0_all_testing_gg = dr0_all_gg[-used_id_gg,]
used_id_chr3_4_gg = sample(dim(dr0_all_testing_gg)[1], 10000)
exp_fitpre_gg = predict(fit_gg, dr0_all_testing_gg[used_id_chr3_4_gg,])
fit_R2pre_gg_svm = R2(exp_fitpre_gg, dr0_all_testing_gg[used_id_chr3_4_gg,55])
exp_fitpre_gg_all_svm = c(exp_fitpre_gg_all_svm, exp_fitpre_gg)
dr0_all_testing_gg_all_svm = c(dr0_all_testing_gg_all_svm, dr0_all_testing_gg[used_id_chr3_4_gg,55])
print(fit_R2pre_gg_svm)
}
fitpre_LM_R2_gg_svm = R2(exp_fitpre_gg_all_svm, dr0_all_testing_gg_all_svm)
fitpre_LM_R2_gg_svm
0.8228507
0.8607901
0.8587032

### All Proximal
exp_fitpre_gg_all_P_svm = c()
dr0_all_testing_gg_all_P_svm = c()
for (gg in c(1:4)){
#gg = 1
dr0_all_gg = dr0_all[label_all_genegroup==gg,]
set.seed(2018)
used_id_gg = sample(dim(dr0_all_gg)[1], 10000)
fit_gg = svm(Gene_Group ~ .-1, data=dr0_all_gg[used_id_gg,c(28:54,55)])
fit_R2_gg = summary(fit_gg)$r.squared
print(fit_R2_gg)
dr0_all_testing_gg = dr0_all_gg[-used_id_gg,]
used_id_chr3_4_gg = sample(dim(dr0_all_testing_gg)[1], 10000)
exp_fitpre_gg = predict(fit_gg, dr0_all_testing_gg[used_id_chr3_4_gg,c(28:54,55)])
fit_R2pre_gg_P_svm = R2(exp_fitpre_gg, dr0_all_testing_gg[used_id_chr3_4_gg,55])
exp_fitpre_gg_all_P_svm = c(exp_fitpre_gg_all_P_svm, exp_fitpre_gg)
dr0_all_testing_gg_all_P_svm = c(dr0_all_testing_gg_all_P_svm, dr0_all_testing_gg[used_id_chr3_4_gg,55])
print(fit_R2pre_gg_P_svm)
}
fitpre_LM_R2_gg_P_svm = R2(exp_fitpre_gg_all_P_svm, dr0_all_testing_gg_all_P_svm)
fitpre_LM_R2_gg_P_svm
0.6368281
0.7177365
0.7206636

### All Distal
exp_fitpre_gg_all_D_svm = c()
dr0_all_testing_gg_all_D_svm = c()
for (gg in c(1:4)){
#gg = 1
dr0_all_gg = dr0_all[label_all_genegroup==gg,]
set.seed(2018)
used_id_gg = sample(dim(dr0_all_gg)[1], 10000)
fit_gg = svm(Gene_Group ~ .-1, data=dr0_all_gg[used_id_gg,c(1:27,55)])
fit_R2_gg = summary(fit_gg)$r.squared
print(fit_R2_gg)
dr0_all_testing_gg = dr0_all_gg[-used_id_gg,]
used_id_chr3_4_gg = sample(dim(dr0_all_testing_gg)[1], 10000)
exp_fitpre_gg = predict(fit_gg, dr0_all_testing_gg[used_id_chr3_4_gg,c(1:27,55)])
fit_R2pre_gg_D_svm = R2(exp_fitpre_gg, dr0_all_testing_gg[used_id_chr3_4_gg,55])
exp_fitpre_gg_all_D_svm = c(exp_fitpre_gg_all_D_svm, exp_fitpre_gg)
dr0_all_testing_gg_all_D_svm = c(dr0_all_testing_gg_all_D_svm, dr0_all_testing_gg[used_id_chr3_4_gg,55])
print(fit_R2pre_gg_D_svm)
}
fitpre_LM_R2_gg_D_svm = R2(exp_fitpre_gg_all_D_svm, dr0_all_testing_gg_all_D_svm)
fitpre_LM_R2_gg_D_svm
0.8045486
0.8499905
0.8451157


png('test.GG.split.SVM.png', width=2000, height=500)
par(mfrow=c(1,4)) 
used_id_gg_plot = c(1:10000)
plot(exp_fitpre_gg_all_svm[used_id_gg_plot], dr0_all_testing_gg_all_svm[used_id_gg_plot], xlim = c(-10, 15), ylim = c(-10, 15), main=round(R2(exp_fitpre_gg_all_svm[used_id_gg_plot], dr0_all_testing_gg_all_svm[used_id_gg_plot]),3))
abline(0,1, col='red')
used_id_gg_plot = c(10001:20000)
plot(exp_fitpre_gg_all_svm[used_id_gg_plot], dr0_all_testing_gg_all_svm[used_id_gg_plot], xlim = c(-10, 15), ylim = c(-10, 15), main=round(R2(exp_fitpre_gg_all_svm[used_id_gg_plot], dr0_all_testing_gg_all_svm[used_id_gg_plot]),3))
abline(0,1, col='red')
used_id_gg_plot = c(20001:30000)
plot(exp_fitpre_gg_all_svm[used_id_gg_plot], dr0_all_testing_gg_all_svm[used_id_gg_plot], xlim = c(-10, 15), ylim = c(-10, 15), main=round(R2(exp_fitpre_gg_all_svm[used_id_gg_plot], dr0_all_testing_gg_all_svm[used_id_gg_plot]),3))
abline(0,1, col='red')
used_id_gg_plot = c(30001:40000)
plot(exp_fitpre_gg_all_svm[used_id_gg_plot], dr0_all_testing_gg_all_svm[used_id_gg_plot], xlim = c(-10, 15), ylim = c(-10, 15), main=round(R2(exp_fitpre_gg_all_svm[used_id_gg_plot], dr0_all_testing_gg_all_svm[used_id_gg_plot]),3))
abline(0,1, col='red')
dev.off()


png('test.GG.split.D.SVM.png', width=2000, height=500)
par(mfrow=c(1,4)) 
used_id_gg_plot = c(1:10000)
plot(exp_fitpre_gg_all_D_svm[used_id_gg_plot], dr0_all_testing_gg_all_D_svm[used_id_gg_plot], xlim = c(-10, 15), ylim = c(-10, 15), main=round(R2(exp_fitpre_gg_all_D_svm[used_id_gg_plot], dr0_all_testing_gg_all_D_svm[used_id_gg_plot]),3))
abline(0,1, col='red')
used_id_gg_plot = c(10001:20000)
plot(exp_fitpre_gg_all_D_svm[used_id_gg_plot], dr0_all_testing_gg_all_D_svm[used_id_gg_plot], xlim = c(-10, 15), ylim = c(-10, 15), main=round(R2(exp_fitpre_gg_all_D_svm[used_id_gg_plot], dr0_all_testing_gg_all_D_svm[used_id_gg_plot]),3))
abline(0,1, col='red')
used_id_gg_plot = c(20001:30000)
plot(exp_fitpre_gg_all_D_svm[used_id_gg_plot], dr0_all_testing_gg_all_D_svm[used_id_gg_plot], xlim = c(-10, 15), ylim = c(-10, 15), main=round(R2(exp_fitpre_gg_all_D_svm[used_id_gg_plot], dr0_all_testing_gg_all_D_svm[used_id_gg_plot]),3))
abline(0,1, col='red')
used_id_gg_plot = c(30001:40000)
plot(exp_fitpre_gg_all_D_svm[used_id_gg_plot], dr0_all_testing_gg_all_D_svm[used_id_gg_plot], xlim = c(-10, 15), ylim = c(-10, 15), main=round(R2(exp_fitpre_gg_all_D_svm[used_id_gg_plot], dr0_all_testing_gg_all_D_svm[used_id_gg_plot]),3))
abline(0,1, col='red')
dev.off()


png('test.GG.split.P.SVM.png', width=2000, height=500)
par(mfrow=c(1,4)) 
used_id_gg_plot = c(1:10000)
plot(exp_fitpre_gg_all_P_svm[used_id_gg_plot], dr0_all_testing_gg_all_P_svm[used_id_gg_plot], xlim = c(-10, 15), ylim = c(-10, 15), main=round(R2(exp_fitpre_gg_all_P_svm[used_id_gg_plot], dr0_all_testing_gg_all_P_svm[used_id_gg_plot]),3))
abline(0,1, col='red')
used_id_gg_plot = c(10001:20000)
plot(exp_fitpre_gg_all_P_svm[used_id_gg_plot], dr0_all_testing_gg_all_P_svm[used_id_gg_plot], xlim = c(-10, 15), ylim = c(-10, 15), main=round(R2(exp_fitpre_gg_all_P_svm[used_id_gg_plot], dr0_all_testing_gg_all_P_svm[used_id_gg_plot]),3))
abline(0,1, col='red')
used_id_gg_plot = c(20001:30000)
plot(exp_fitpre_gg_all_P_svm[used_id_gg_plot], dr0_all_testing_gg_all_P_svm[used_id_gg_plot], xlim = c(-10, 15), ylim = c(-10, 15), main=round(R2(exp_fitpre_gg_all_P_svm[used_id_gg_plot], dr0_all_testing_gg_all_P_svm[used_id_gg_plot]),3))
abline(0,1, col='red')
used_id_gg_plot = c(30001:40000)
plot(exp_fitpre_gg_all_P_svm[used_id_gg_plot], dr0_all_testing_gg_all_P_svm[used_id_gg_plot], xlim = c(-10, 15), ylim = c(-10, 15), main=round(R2(exp_fitpre_gg_all_P_svm[used_id_gg_plot], dr0_all_testing_gg_all_P_svm[used_id_gg_plot]),3))
abline(0,1, col='red')
dev.off()

png('test.GG.SVM.png', width=1500, height=500)
par(mfrow=c(1,3)) 
plot(exp_fitpre_gg_all_svm, dr0_all_testing_gg_all_svm, xlim = c(-10, 15), ylim = c(-10, 15), main=round(R2(exp_fitpre_gg_all_svm, dr0_all_testing_gg_all_svm),3))
abline(0,1, col='red')
plot(exp_fitpre_gg_all_P_svm, dr0_all_testing_gg_all_P_svm, xlim = c(-10, 15), ylim = c(-10, 15), main=round(R2(exp_fitpre_gg_all_P_svm, dr0_all_testing_gg_all_P_svm),3))
abline(0,1, col='red')
plot(exp_fitpre_gg_all_D_svm, dr0_all_testing_gg_all_D_svm, xlim = c(-10, 15), ylim = c(-10, 15), main=round(R2(exp_fitpre_gg_all_D_svm, dr0_all_testing_gg_all_D_svm),3))
abline(0,1, col='red')
dev.off()










### All RF
modelRF = randomForest(dr0_all[used_id_chr1_2,-55], dr0_all[used_id_chr1_2,55], ntree = 500, mtry = 6, importance = TRUE)
fit_RF_pre = predict(modelRF, dr0_all_testing[used_id_chr3_4,-55])
fit_RF_R2pre = R2(fit_RF_pre, dr0_all_testing[used_id_chr3_4,55])
fit_RF_R2pre
0.6164033
### Proximal RF
modelRF_P = randomForest(dr0_all[used_id_chr1_2,c(28:54)], dr0_all[used_id_chr1_2,55], ntree = 500, mtry = 6, importance = TRUE)
fit_RF_P_pre = predict(modelRF_P, dr0_all_testing[used_id_chr3_4,c(28:54)])
fit_RF_P_R2pre = R2(fit_RF_P_pre, dr0_all_testing[used_id_chr3_4,55])
fit_RF_P_R2pre
0.4895566
### Distal RF
modelRF_D = randomForest(dr0_all[used_id_chr1_2,c(1:27)], dr0_all[used_id_chr1_2,55], ntree = 500, mtry = 6, importance = TRUE)
fit_RF_D_pre = predict(modelRF_D, dr0_all_testing[used_id_chr3_4,-55])
fit_RF_D_R2pre = R2(fit_RF_D_pre, dr0_all_testing[used_id_chr3_4,55])
fit_RF_D_R2pre
0.4747674

png('test.RF.png', width=1500, height=500)
par(mfrow=c(1,3)) 
plot(fit_RF_pre, dr0_all_testing[used_id_chr3_4,55], xlim = c(-10, 15), ylim = c(-10, 15))
abline(0,1, col='red')
plot(fit_RF_P_pre, dr0_all_testing[used_id_chr3_4,55], xlim = c(-10, 15), ylim = c(-10, 15))
abline(0,1, col='red')
plot(fit_RF_D_pre, dr0_all_testing[used_id_chr3_4,55], xlim = c(-10, 15), ylim = c(-10, 15))
abline(0,1, col='red')
dev.off()




### All NN
f = 'Gene_Group ~ '
for (cn in colnames_state[1:53]){
	f = paste(f, cn, '+', sep=' ')
}
f = paste(f, colnames_state[54], sep=' ')

used_id_chr1_2_nn = sample(dim(dr0_all)[1], 1000)
nn = neuralnet(f, data=dr0_all[used_id_chr1_2_nn,], hidden=c(20,10),linear.output=T)
pr.nn = compute(nn, dr0_all_testing[used_id_chr3_4,-55])
fitpre_NN_R2 = R2(pr.nn$net.result, dr0_all_testing[used_id_chr3_4,55])
fitpre_NN_R2





fit_lm = lm(dr0_all[,55]~dr0_all[,-55]-1)
eRP = fit_lm$coefficients
eRP_mat_D_all = eRP[1:27]
eRP_mat_P_all = eRP[28:54]
eRP_mat = cbind(eRP_mat_D_all, eRP_mat_P_all, eRP_mat)
dr_kmeans_plot_cor0 = cor(dr0_all, method='pearson')[-55,55]
dr_kmeans_plot_cor_all = cbind(dr_kmeans_plot_cor0[1:27], dr_kmeans_plot_cor0[28:54], dr_kmeans_plot_cor_all[1:27,], dr_kmeans_plot_cor_all[28:54,])


corlo_lim_scale = max(c(max(abs(dr_kmeans_plot_cor_all)), max(abs(dr_kmeans_plot_cor_all))))
breaksList_scale = seq(-corlo_lim_scale, corlo_lim_scale, by = 0.01)
print(breaksList_scale)
my_colorbar_scale=colorRampPalette(c('blue', 'white', 'red'))(n = length(breaksList_scale))

library(pheatmap)
nr=20
colnames(dr_kmeans_plot_cor_all) = c('D0', 'P0', 'D1', 'D2', 'D3', 'D4', 'P1', 'P2', 'P3', 'P4')
rownames(dr_kmeans_plot_cor_all) = c(0:26)
png(paste('kmean.', toString(nr), '.pcor.png', sep=''), width=300)
pheatmap(dr_kmeans_plot_cor_all[rev(fit$order),-c(1,2)], color=my_colorbar_scale, breaks = breaksList_scale, cluster_rows = F, cluster_cols = F,annotation_names_row=FALSE,annotation_names_col=FALSE,show_rownames=T,show_colnames=T)
dev.off()
png(paste('kmean.', toString(nr), '.pcor0.png', sep=''), width=150)
pheatmap(dr_kmeans_plot_cor_all[rev(fit$order),1:2], color=my_colorbar_scale, breaks = breaksList_scale, cluster_rows = F, cluster_cols = F,annotation_names_row=FALSE,annotation_names_col=FALSE,show_rownames=T,show_colnames=T)
dev.off()
png(paste('kmean.', toString(nr), '.pcor0.cater.png', sep=''), width=150)
pheatmap(cater_mat[rev(fit$order),], color=my_colorbar_scale, breaks = breaksList_scale, cluster_rows = F, cluster_cols = F,annotation_names_row=FALSE,annotation_names_col=FALSE,show_rownames=T,show_colnames=T)
dev.off()
png(paste('kmean.', toString(nr), '.pcor.dif.pdf', sep=''), width=200)
pheatmap(dr_kmeans_plot_cor_all[rev(fit$order),c(3:6)]-dr_kmeans_plot_cor_all[rev(fit$order),c(7:10)], color=my_colorbar_scale, breaks = breaksList_scale, cluster_rows = F, cluster_cols = F,annotation_names_row=FALSE,annotation_names_col=FALSE,show_rownames=T,show_colnames=T)
dev.off()
colnames(eRP_mat) = c('D0', 'P0', 'D1', 'D2', 'D3', 'D4', 'P1', 'P2', 'P3', 'P4')
rownames(eRP_mat) = c(0:26)
min(eRP_mat[!is.na(eRP_mat)])
png(paste('kmean.', toString(nr), '.eRP.png', sep=''), width=300)
pheatmap(eRP_mat[rev(fit$order),-c(1,2)], color=my_colorbar_scale, cluster_rows = F, cluster_cols = F,annotation_names_row=FALSE,annotation_names_col=FALSE,show_rownames=T,show_colnames=T)
dev.off()
png(paste('kmean.', toString(nr), '.eRP0.png', sep=''), width=150)
pheatmap(eRP_mat[rev(fit$order),c(1,2)], color=my_colorbar_scale, cluster_rows = F, cluster_cols = F,annotation_names_row=FALSE,annotation_names_col=FALSE,show_rownames=T,show_colnames=T)
dev.off()

eRP_mat_center = apply(eRP_mat, 2, function(x) x-x[1])
corlo_lim_scale_eRP = max(abs(eRP_mat_center[!is.na(eRP_mat)]))
breaksList_scale_eRP = seq(-corlo_lim_scale_eRP, corlo_lim_scale_eRP, by = 1)
my_colorbar_scale_eRP =colorRampPalette(c('blue', 'white', 'red'))(n = length(breaksList_scale_eRP))
colnames(eRP_mat_center) = c('D0', 'P0', 'D1', 'D2', 'D3', 'D4', 'P1', 'P2', 'P3', 'P4')

png(paste('kmean.', toString(nr), '.eRPcenter.png', sep=''), width=300)
pheatmap(eRP_mat_center[rev(fit$order),-c(1,2)], color=my_colorbar_scale_eRP, breaks = breaksList_scale_eRP, cluster_rows = F, cluster_cols = F,annotation_names_row=FALSE,annotation_names_col=FALSE,show_rownames=T,show_colnames=T)
dev.off()
png(paste('kmean.', toString(nr), '.eRPcenter0.png', sep=''), width=150)
pheatmap(eRP_mat_center[rev(fit$order),c(1,2)], color=my_colorbar_scale_eRP, breaks = breaksList_scale_eRP, cluster_rows = F, cluster_cols = F,annotation_names_row=FALSE,annotation_names_col=FALSE,show_rownames=T,show_colnames=T)
dev.off()


eRP_mat_center = apply(eRP_mat, 2, function(x) x-x[1])
eRP_mat_center_p = eRP_mat_center
eRP_mat_center_p[eRP_mat_center>0] = log2(1/eRP_mat_center[eRP_mat_center>0])
eRP_mat_center_p[eRP_mat_center<0] = -log2(-1/eRP_mat_center[eRP_mat_center<0])

lim0=1
eRP_mat_center_p[eRP_mat_center_p>lim0]=lim0
eRP_mat_center_p[eRP_mat_center_p<(-lim0)]=-lim0
corlo_lim_scale_eRP = max(abs(eRP_mat_center_p[-1,][!is.na(eRP_mat[-1,])]))
breaksList_scale_eRP = seq(-corlo_lim_scale_eRP, corlo_lim_scale_eRP, by = 0.01)
my_colorbar_scale_eRP =colorRampPalette(c('white', 'red', 'white'))(n = length(breaksList_scale_eRP))
colnames(eRP_mat_center) = c('D0', 'P0', 'D1', 'D2', 'D3', 'D4', 'P1', 'P2', 'P3', 'P4')

png(paste('kmean.', toString(nr), '.eRPcenter_p.png', sep=''), width=300)
pheatmap(eRP_mat_center_p[rev(fit$order),-c(1,2)], color=my_colorbar_scale_eRP, breaks = breaksList_scale_eRP, cluster_rows = F, cluster_cols = F,annotation_names_row=FALSE,annotation_names_col=FALSE,show_rownames=T,show_colnames=T)
dev.off()
png(paste('kmean.', toString(nr), '.eRPcenter0_p.png', sep=''), width=150)
pheatmap(eRP_mat_center_p[rev(fit$order),c(1,2)], color=my_colorbar_scale_eRP, breaks = breaksList_scale_eRP, cluster_rows = F, cluster_cols = F,annotation_names_row=FALSE,annotation_names_col=FALSE,show_rownames=T,show_colnames=T)
dev.off()


pdf(paste('kmean.', toString(nr), '.pcor.pdf', sep=''), width=5)
pheatmap(dr_kmeans_plot_cor_all[rev(fit$order),-c(1,2)], color=my_colorbar_scale, breaks = breaksList_scale, cluster_rows = F, cluster_cols = F,annotation_names_row=FALSE,annotation_names_col=FALSE,show_rownames=T,show_colnames=T)
dev.off()
pdf(paste('kmean.', toString(nr), '.pcor0.pdf', sep=''), width=3)
pheatmap(dr_kmeans_plot_cor_all[rev(fit$order),1:2], color=my_colorbar_scale, breaks = breaksList_scale, cluster_rows = F, cluster_cols = F,annotation_names_row=FALSE,annotation_names_col=FALSE,show_rownames=T,show_colnames=T)
dev.off()
pdf(paste('kmean.', toString(nr), '.pcor0.cater.pdf', sep=''), width=3)
pheatmap(cater_mat[rev(fit$order),], color=my_colorbar_scale, breaks = breaksList_scale, cluster_rows = F, cluster_cols = F,annotation_names_row=FALSE,annotation_names_col=FALSE,show_rownames=T,show_colnames=T)
dev.off()
pdf(paste('kmean.', toString(nr), '.pcor.dif.pdf', sep=''), width=4)
pheatmap(dr_kmeans_plot_cor_all[rev(fit$order),c(3:6)]-dr_kmeans_plot_cor_all[rev(fit$order),c(7:10)], color=my_colorbar_scale, breaks = breaksList_scale, cluster_rows = F, cluster_cols = F,annotation_names_row=FALSE,annotation_names_col=FALSE,show_rownames=T,show_colnames=T)
dev.off()
pdf(paste('kmean.', toString(nr), '.eRP.pdf', sep=''), width=5)
pheatmap(eRP_mat[rev(fit$order),], color=my_colorbar_scale, cluster_rows = F, cluster_cols = F,annotation_names_row=FALSE,annotation_names_col=FALSE,show_rownames=T,show_colnames=T)
dev.off()


write.table(dr_kmeans_plot_cor_all, 'dr_kmeans_plot_cor_all.txt', quote=F, col.names=F, row.names=F, sep='\t')


source('~/group/projects/vision/createGenomeTracks.R')
pdf('pknorm_2_16lim_ref1mo_0424_heatmap.sig_sort.hclust.2kb.pdf', width=4)
createHeatmap('~/group/projects/vision/pknorm_2_16lim_ref1mo_sample1mo_0424/pknorm_2_16lim_ref1mo_0424_lesshet.para0', sortstate=TRUE)
dev.off()







