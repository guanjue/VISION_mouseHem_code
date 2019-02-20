for i in {1..12}
do
#echo $i
#cat 'vision_rna.chr'*'.'*'.'$i'.gene_ccRE.txt' > 'vision_rna.gene_ccRE.c'$i'.txt'
#cat 'vision_rna.gene_ccRE.c'$i'.txt' | awk -F '\t' -v OFS='\t' '{if ($8==1) print $0}' > 'vision_rna.gene_ccRE.c'$i'.selected.txt'
cat 'vision_rna.gene_ccRE.c'$i'.selected.txt' | awk -F '\t' -v OFS='\t' '{if ($3=="protein_coding") print $0}' > 'vision_rna.gene_ccRE.c'$i'.selected.protein_coding.txt'
cat 'vision_rna.gene_ccRE.c'$i'.selected.txt' | awk -F '\t' -v OFS='\t' '{if ($3!="protein_coding") print $0}' > 'vision_rna.gene_ccRE.c'$i'.selected.other.txt'
wc -l 'vision_rna.gene_ccRE.c'$i'.selected.txt'
done


R

d = read.table('vision_rna.gene_ccRE.c1.selected.protein_coding.txt', header=F)
tss_vec = apply(d,1,function(x) paste(x[1],x[2],collapes='_'))


d=read.table('vision_rna.gene_ccRE.selected.protein_coding.reprod_count.txt', header=F)
d_reprod = d[d[,9]==12,]

tss_d_reprod_vec = apply(d_reprod,1,function(x) paste(x[1],x[2],collapes='_'))
pdf('hist_reprod.pdf')
hist(table(tss_d_reprod_vec), breaks=50)
dev.off()



tss_adjr2 = c()
distal_adjr2 = c()
tss_distal_adjr2 = c()

k = 0

for (h in c(1:19)){
for (j in c(1:4)){
for (i in c(1:12)){
	k = k+1
	filename = paste('vision_rna.chr', toString(h), '.', toString(j), '.', toString(i), '.Rdat', sep='')
	load(filename)
	fit = lm(rt$y~log2(rt$x+0.001))
	fit_tss = lm(rt$y~log2(rt$x0+0.001))
	fit_tss_distal = lm(rt$y~log2(rt$x+0.001)+log2(rt$x0+0.001))
	tss_adjr2[k] = summary(fit_tss)$adj.r.squared
	distal_adjr2[k] = summary(fit)$adj.r.squared
	tss_distal_adjr2[k] = summary(fit_tss_distal)$adj.r.squared
	#print(summary(fit_tss)$adj.r.squared)
	#print(summary(fit)$adj.r.squared)
}
}
print(h)
}


pdf('tss_vs_distal.pdf')
boxplot(cbind(tss_adjr2, distal_adjr2, tss_distal_adjr2))
dev.off()

