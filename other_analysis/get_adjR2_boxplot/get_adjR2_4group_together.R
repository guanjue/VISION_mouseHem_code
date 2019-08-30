setwd('/gpfs/group/yzz2/default/legacy/group/projects/vision/rna')



x0_mat = c()
x_mat = c()
y_mat = c()

tss_adjr2 = c()
distal_adjr2 = c()
tss_distal_adjr2 = c()
k=0

for (h in c(1:19, 'X')){
print(h)
k = k+1
for (j in 1:4){

	i=1
	filename = paste('vision_rna_tss2k_ccreunit_folder_used0509/vision_rna_tss2k_ccreunit.chr', toString(h), '.', toString(j), '.', toString(i), '.Rdat', sep='')
	load(filename)

	y_mat = c(y_mat, rt$y)
	x_mat = rbind(x_mat, rt$x)
	x0_mat = rbind(x0_mat, rt$x0)
	#print(summary(fit_tss)$adj.r.squared)
	#print(summary(fit)$adj.r.squared)

}

fit = lm(y_mat~log2(x_mat+0.001))
fit_tss = lm(y_mat~log2(x0_mat+0.001))
fit_tss_distal = lm(y_mat~log2(x_mat+0.001)+log2(x0_mat+0.001))

tss_adjr2[k] = summary(fit_tss)$adj.r.squared
distal_adjr2[k] = summary(fit)$adj.r.squared
tss_distal_adjr2[k] = summary(fit_tss_distal)$adj.r.squared

}


pdf(paste('tss_vs_distal.all4group.pdf', sep=''), width=3.5, height=4)
boxplot(cbind(tss_adjr2, distal_adjr2, tss_distal_adjr2), col=c('red', 'blue', 'orange'), names = c('P', 'D', 'P & D'), ylim=c(0,1))
dev.off()


