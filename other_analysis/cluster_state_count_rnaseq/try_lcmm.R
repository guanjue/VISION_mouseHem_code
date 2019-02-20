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
chr_num = c(1,3, 'X')
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
used_id_chr1_2 =  sample(dim(dr0_all)[1], 1000)
used_id_chr3_4 =  sample(dim(dr0_all_testing)[1], 40000)
used_id_chr1_2_nn =  sample(dim(dr0_all)[1], 1000)



fit = lm(Gene_Group ~ .-1, data=dr0_all[used_id_chr1_2,])

summary(fit)

d2 <-lcmm(Gene_Group ~ P1+P2,mixture=~P1+P2,subject='P1',ng=2,idiag=TRUE,link="linear", data=dr0_all[used_id_chr1_2,])



summary(d2) 



