setwd('/gpfs/group/yzz2/default/legacy/group/projects/vision/rna')

for (gg in 1:4){
tss_adjr2 = c()
distal_adjr2 = c()
tss_distal_adjr2 = c()

k = 0

for (h in c(1:19)){
for (j in c(gg)){
#for (i in c(1:12)){
	i=1
	k = k+1
	filename = paste('vision_rna_tss2k_ccreunit_folder_used0509/vision_rna_tss2k_ccreunit.chr', toString(h), '.', toString(j), '.', toString(i), '.Rdat', sep='')
	load(filename)
	fit = lm(rt$y~log2(rt$x+0.001))
	fit_tss = lm(rt$y~log2(rt$x0+0.001))
	fit_tss_distal = lm(rt$y~log2(rt$x+0.001)+log2(rt$x0+0.001))
	tss_adjr2[k] = summary(fit_tss)$adj.r.squared
	distal_adjr2[k] = summary(fit)$adj.r.squared
	tss_distal_adjr2[k] = summary(fit_tss_distal)$adj.r.squared
	#print(summary(fit_tss)$adj.r.squared)
	#print(summary(fit)$adj.r.squared)
#}
}
print(h)
}

pdf(paste('tss_vs_distal.', gg, '.pdf', sep=''), width=3.5, height=4)
boxplot(cbind(tss_adjr2, distal_adjr2, tss_distal_adjr2), col=c('red', 'blue', 'orange'), names = c('P', 'D', 'P & D'), ylim=c(0,1))
dev.off()
print(gg)
}

