setwd('/gpfs/group/yzz2/default/legacy/group/projects/vision/rna')


adjR2_fun = function(x, xp, k) {
	R2 = 1-sum((x-xp)^2)/sum((x-mean(x))^2)
	n = length(x)
	R2adj = 1 - (1-R2)*(n-1)/(n-(k+1))
	return(R2adj)
}

R2_fun = function(x, xp, k) {
	R2 = 1-sum((x-xp)^2)/sum((x-mean(x))^2)
	return(R2)
}


for (gg in 1:4){
adjR2_vec_all = c()
adjR2_vec_tss = c()
adjR2_vec_dist = c()

k = 0
smallnum = 0.00

y_g = c()
x_g = c()
x0_g = c()
for (h in c(1:19, 'X')){
	for (j in c(gg)){
	i=1
		k = k+1
		filename = paste('vision_rna_tss2k_ccreunit_folder_used0509/vision_rna_tss2k_ccreunit.chr', toString(h), '.', toString(j), '.', toString(i), '.Rdat', sep='')
		load(filename)
		y_all = rt$y
		x_all = rt$nx
		x0_all = rt$x0
		one_ct_num = length(y_all)/12
		y_chrtmp = c()
		x_chrtmp = c()
		x0_chrtmp = c()
		### get all info in chr
		for (l in c(1:12)){
			cttmp_id = c(1:one_ct_num) + (l-1)*one_ct_num
			y_cttmp = y_all[cttmp_id]
			x_cttmp = (x_all[cttmp_id,]+smallnum)
			x0_cttmp = (x0_all[cttmp_id,]+smallnum)	
			y_chrtmp = cbind(y_chrtmp, y_cttmp)
			x_chrtmp = cbind(x_chrtmp, x_cttmp)
			x0_chrtmp = cbind(x0_chrtmp, x0_cttmp)
		}
		### rbind
		y_g = rbind(y_g, y_chrtmp)
		x_g = rbind(x_g, x_chrtmp)
		x0_g = rbind(x0_g, x0_chrtmp)
	}
	print(h)
}

### get performance
for (l in c(1:12)){
	print(l)
	### train vs test y
	y_g_testing = y_g[,l]
	y_g_training = as.vector(y_g[,-l])
	### test x & x0
	used_id = c(1:27) + (l-1)*27
	x_g_testing = x_g[,used_id]
	x0_g_testing = x0_g[,used_id]
	### train x & x0
	x_g_training_all = x_g[,-used_id]
	x_g_training = c()
	x0_g_training_all = x0_g[,-used_id]
	x0_g_training = c()
	for (m in c(1:11)){
		used_id_m = c(1:27) + (m-1)*27
		x_g_training = rbind(x_g_training, x_g_training_all[,used_id_m])
		x0_g_training = rbind(x0_g_training, x0_g_training_all[,used_id_m])
	}
	### train model
	data_train_all = as.data.frame(cbind(y_g_training, x_g_training, x0_g_training))
	colnames(data_train_all) = c('y',c(1:(2*dim(x_g_training)[2])))
	fit_all = lm(y~., data=data_train_all)
	data_train_tss = as.data.frame(cbind(y_g_training, x0_g_training))
	colnames(data_train_tss) = c('y',c(1:dim(x_g_training)[2]))
	fit_tss = lm(y~., data=data_train_tss)
	data_train_dist = as.data.frame(cbind(y_g_training, x_g_training))
	colnames(data_train_dist) = c('y',c(1:dim(x_g_training)[2]))
	fit_dist = lm(y~., data=data_train_dist)
	### testing
	data_test_all = as.data.frame(cbind(x_g_testing, x0_g_testing))
	colnames(data_test_all) = c(c(1:(2*dim(x_g_training)[2])))
	y_test_pred_all = predict(fit_all, newdata=data_test_all)
	adjR2_all = adjR2_fun(y_g_testing, y_test_pred_all, 2*dim(x_g_training)[2])
	data_test_tss = as.data.frame(cbind(x0_g_testing))
	colnames(data_test_tss) = c(c(1:dim(x_g_training)[2]))
	y_test_pred_tss = predict(fit_tss, newdata=data_test_tss)
	adjR2_tss = adjR2_fun(y_g_testing, y_test_pred_tss, dim(x_g_training)[2])
	data_test_dist = as.data.frame(cbind(x_g_testing))
	colnames(data_test_dist) = c(c(1:dim(x_g_training)[2]))
	y_test_pred_dist = predict(fit_dist, newdata=data_test_dist)
	adjR2_dist = adjR2_fun(y_g_testing, y_test_pred_dist, dim(x_g_training)[2])
	###
	adjR2_vec_all[l] = adjR2_all
	adjR2_vec_tss[l] = adjR2_tss
	adjR2_vec_dist[l] = adjR2_dist
}

pdf(paste('tss_vs_distal.', gg, '.pdf', sep=''), width=3, height=3.5)
boxplot(cbind(adjR2_vec_tss, adjR2_vec_dist, adjR2_vec_all), col=c('red', 'blue', 'orange'), names = c('P', 'D', 'P & D'), ylim=c(-0.2,1))
dev.off()
print(gg)
}

