setwd('/gpfs/group/yzz2/default/legacy/group/projects/vision/rna')


adjR2_fun = function(x, xp, k) {
	R2 = 1-sum((x-xp)^2)/sum((x-mean(x))^2)
	n = length(x)
	R2adj = 1 - (1-R2)*(n-1)/(n-(k+1))
	return(R2adj)
}


for (gg in 1:4){
tss_adjr2 = c()
distal_adjr2 = c()
tss_distal_adjr2 = c()

k = 0
smallnum = 0.001

for (h in c(1:19, 'X')){
for (j in c(gg)){
#for (i in c(1)){
	k = k+1
	filename = paste('vision_rna_tss2k_ccreunit_folder_used0509/vision_rna_tss2k_ccreunit.chr', toString(h), '.', toString(j), '.', toString(1), '.Rdat', sep='')
	load(filename)
	y_all = rt$y
	x_all = rt$x
	x0_all = rt$x0
	one_ct_num = length(y_all)/12
for (i in c(1:12)){
	test_id = c(1:one_ct_num) + (i-1)*one_ct_num
	y_train = y_all[-test_id]
	x_train = log2(x_all[-test_id,]+smallnum)
	x0_train = log2(x0_all[-test_id,]+smallnum)
	y_test = y_all[test_id]
	x_test = log2(x_all[test_id,]+smallnum)
	x0_test = log2(x0_all[test_id,]+smallnum)
	### training
	training_data = as.data.frame(cbind(y_train, x_train))
	colnames(training_data) = c('y', 1:dim(x_train)[2])
	fit = lm(y~., data=training_data)
	training_data_tss = as.data.frame(cbind(y_train, x0_train))
	colnames(training_data_tss) = c('y', 1:dim(x0_train)[2])
	fit_tss = lm(y~., data=training_data_tss)
	training_data_tss_distal = as.data.frame(cbind(y_train, x_train, x0_train))
	colnames(training_data_tss_distal) = c('y', 1:(dim(x_train)[2]+dim(x0_train)[2]))
	fit_tss_distal = lm(y~., data=training_data_tss_distal)
	### testing
	testing_data = as.data.frame(cbind(x_test))
	colnames(testing_data) = c(1:dim(testing_data)[2])
	y_test_pred = predict(fit, newdata=testing_data)
	adjR2 = adjR2_fun(y_test, y_test_pred, dim(testing_data)[2])

	testing_data_tss = as.data.frame(cbind(x0_test))
	colnames(testing_data_tss) = c(1:dim(testing_data_tss)[2])
	y_test_pred_tss = predict(fit_tss, newdata=testing_data_tss)
	adjR2_tss = adjR2_fun(y_test, y_test_pred_tss, dim(testing_data_tss)[2])

	testing_data_tss_dist = as.data.frame(cbind(x_test, x0_test))
	colnames(testing_data_tss_dist) = c(1:dim(testing_data_tss_dist)[2])
	y_test_pred_tss_dist = predict(fit, newdata=testing_data_tss_dist)
	adjR2_tss_dist = adjR2_fun(y_test, y_test_pred_tss_dist, dim(testing_data_tss_dist)[2])

	tss_adjr2[k] = adjR2_tss
	distal_adjr2[k] = adjR2
	tss_distal_adjr2[k] = adjR2_tss_dist
	#print(summary(fit_tss)$adj.r.squared)
	#print(summary(fit)$adj.r.squared)
}
}
print(h)
}

pdf(paste('tss_vs_distal.', gg, '.pdf', sep=''), width=3.5, height=4)
boxplot(cbind(tss_adjr2, distal_adjr2, tss_distal_adjr2), col=c('red', 'blue', 'orange'), names = c('P', 'D', 'P & D'), ylim=c(-1,1))
dev.off()
print(gg)
}

