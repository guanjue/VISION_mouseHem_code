library(randomForest)
library(pheatmap)
library(e1071)
library(neuralnet)

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
chr_num = c(1:19, 'X')
chr_num = c(1:2)
for (i in c(1:4)){
	print(i)
	dr0 = c()
	for (j in c(1:12)){
		print(j)
		for (k in chr_num){
		print(k)
		load(paste('../vision_rna_tss2k_ccreunit.chr', k, '.', i, '.', j, '.Rdat', sep=''))
		dr0 = rbind(dr0, cbind(rt$x, rt$x0, rt$y))
		}
	}
	dr_kmeans_plot_cor = cor(dr0, method='pearson')
	dr_kmeans_plot_cor_all = cbind(dr_kmeans_plot_cor_all, dr_kmeans_plot_cor[-55,55])
	dr0_all = rbind(dr0_all, dr0)
}



dr0_all_testing = c()
chr_num = c(1:19, 'X')
chr_num_testing = c(3:4)
for (i in c(1:4)){
	print(i)
	dr0 = c()
	for (j in c(1:12)){
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
used_id_chr1_2 =  sample(dim(dr0_all)[1], 10000)
used_id_chr3_4 =  sample(dim(dr0_all_testing)[1], 10000)

### All
fit = lm(Gene_Group ~ .-1, data=dr0_all[used_id_chr1_2,])
fit_R2 = summary(fit)$r.squared
fit_R2
exp_fitpre = predict(fit, as.data.frame(dr0_all_testing[used_id_chr3_4,-55]))
fit_R2pre = R2(exp_fitpre, dr0_all_testing[used_id_chr3_4,55])
fit_R2pre
0.5218587
### Proximal
fit_P = lm(Gene_Group ~ .-1, data=dr0_all[used_id_chr1_2,28:55])
fit_P_R2 = summary(fit_P)$r.squared
fit_P_R2
exp_fitpre_P = predict(fit_P, as.data.frame(dr0_all_testing[used_id_chr3_4,28:54]))
fit_R2pre_P = R2(exp_fitpre_P, dr0_all_testing[used_id_chr3_4,55])
fit_R2pre_P
0.4716412
### Distal
fit_D = lm(Gene_Group ~ .-1, data=dr0_all[used_id_chr1_2,c(1:27,55)])
fit_D_R2 = summary(fit_D)$r.squared
fit_D_R2
exp_fitpre_D = predict(fit_D, as.data.frame(dr0_all_testing[used_id_chr3_4,c(1:27,55)]))
fit_R2pre_D = R2(exp_fitpre_D, dr0_all_testing[used_id_chr3_4,55])
fit_R2pre_D
0.2822599

png('test.LM.png', width=1500, height=500)
par(mfrow=c(1,3)) 
plot(exp_fitpre, dr0_all_testing[used_id_chr3_4,55], xlim = c(-10, 15), ylim = c(-10, 15))
abline(0,1, col='red')
plot(exp_fitpre_P, dr0_all_testing[used_id_chr3_4,55], xlim = c(-10, 15), ylim = c(-10, 15))
abline(0,1, col='red')
plot(exp_fitpre_D, dr0_all_testing[used_id_chr3_4,55], xlim = c(-10, 15), ylim = c(-10, 15))
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


### All SVR
modelsvm = svm(Gene_Group ~ ., data=dr0_all[used_id_chr1_2,])
fitpre_svr = predict(modelsvm, as.data.frame(dr0_all_testing[used_id_chr3_4,-55]))
fitpre_svrR2 = R2(fitpre_svr, dr0_all_testing[used_id_chr3_4,55])
fitpre_svrR2
0.6016268
### Promixal SVR
modelsvm_P = svm(Gene_Group ~ ., data=dr0_all[used_id_chr1_2,28:55])
fitpre_P_svr = predict(modelsvm_P, as.data.frame(dr0_all_testing[used_id_chr3_4,28:55]))
fitpre_P_svrR2 = R2(fitpre_P_svr, dr0_all_testing[used_id_chr3_4,55])
fitpre_P_svrR2
0.4352918
### Distal SVR
modelsvm_D = svm(Gene_Group ~ ., data=dr0_all[used_id_chr1_2,c(1:27,55)])
fitpre_D_svr = predict(modelsvm_D, as.data.frame(dr0_all_testing[used_id_chr3_4,c(1:27,55)]))
fitpre_D_svrR2 = R2(fitpre_D_svr, dr0_all_testing[used_id_chr3_4,55])
fitpre_P_svrR2
0.4535795

png('test.SVR.png', width=1500, height=500)
par(mfrow=c(1,3)) 
plot(fitpre_svr, dr0_all_testing[used_id_chr3_4,55], xlim = c(-10, 15), ylim = c(-10, 15))
abline(0,1, col='red')
plot(fitpre_P_svr, dr0_all_testing[used_id_chr3_4,55], xlim = c(-10, 15), ylim = c(-10, 15))
abline(0,1, col='red')
plot(fitpre_D_svr, dr0_all_testing[used_id_chr3_4,55], xlim = c(-10, 15), ylim = c(-10, 15))
abline(0,1, col='red')
dev.off()


### All NN
f = 'Gene_Group ~ '
for (cn in colnames_state[1:53]){
	f = paste(f, cn, '+', sep=' ')
}
f = paste(f, colnames_state[54], sep=' ')

used_id_chr1_2_nn = sample(dim(dr0_all)[1], 1000)
nn = neuralnet(f, data=dr0_all[used_id_chr1_2_nn,], hidden=c(20,5),linear.output=T)
pr.nn = compute(nn, dr0_all_testing[used_id_chr3_4,-55])
fitpre_NN_R2 = R2(pr.nn$net.result, dr0_all_testing[used_id_chr3_4,55])
fitpre_NN_R2




fit_c = lm(cater_sig_label[,55]~cater_sig_label[,-55]-1)
adjR2_c_P = 0.9087
fit_c_P = lm(cater_sig_label[,55]~cater_sig_label[,28:54]-1)
adjR2_c_P = 0.9039
fit_c_D = lm(cater_sig_label[,55]~cater_sig_label[,1:27]-1)
adjR2_c_D = 0.8533



library(randomForest)
set.seed(2018)
used_id = sample(dim(cater_sig_label)[1],10000)
cater_sig_label_df = as.data.frame(cater_sig_label)
cater_sig_label_train = cater_sig_label_df[used_id,]
cater_sig_label_left = cater_sig_label_df[-used_id,]
used_id_test = sample(dim(cater_sig_label_left)[1],10000)
cater_sig_label_test = cater_sig_label_left[used_id_test,]

modelRF = randomForest(cater_sig_label_train[,-55], as.factor(cater_sig_label_train[,55]), ntree = 500, mtry = 6, importance = TRUE)
predTrain = predict(modelRF, cater_sig_label_test[,-55], type = "class")
# Checking classification accuracy
table(predTrain, cater_sig_label_test$Gene_Group)  
accuracy = mean(predTrain == cater_sig_label_test$Gene_Group) 
0.6781

used_col = c(1:27)
modelRF_D = randomForest(cater_sig_label_train[,used_col], as.factor(cater_sig_label_train[,55]), ntree = 500, mtry = 6, importance = TRUE)
predTrain = predict(modelRF_D, cater_sig_label_test[,used_col], type = "class")
# Checking classification accuracy
table(predTrain, cater_sig_label_test$Gene_Group)  
accuracy_D = mean(predTrain == cater_sig_label_test$Gene_Group) 
0.6278

used_col = c(28:54)
modelRF_P = randomForest(cater_sig_label_train[,used_col], as.factor(cater_sig_label_train[,55]), ntree = 500, mtry = 6, importance = TRUE)
predTrain = predict(modelRF_P, cater_sig_label_test[,used_col], type = "class")
# Checking classification accuracy
table(predTrain, cater_sig_label_test$Gene_Group)  
accuracy_P = mean(predTrain == cater_sig_label_test$Gene_Group) 
0.5303

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


write.table(dr_kmeans_plot_cor_all, 'dr_kmeans_plot_cor_all.txt', quote=F, col.names=F, row.names=F, sep='\t')


source('~/group/projects/vision/createGenomeTracks.R')
pdf('pknorm_2_16lim_ref1mo_0424_heatmap.sig_sort.hclust.2kb.pdf', width=4)
createHeatmap('~/group/projects/vision/pknorm_2_16lim_ref1mo_sample1mo_0424/pknorm_2_16lim_ref1mo_0424_lesshet.para0', sortstate=TRUE)
dev.off()











