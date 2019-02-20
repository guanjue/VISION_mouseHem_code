

rna = (read.table('rnaTPM.txt', header=T))
ccREs = (read.table('vision_cres.txt', header=F, sep='\t'))


library(LSD)
png('sd_vs_mean.png')
mean_vec = apply(rna[,5:16], 1, mean)
sd_vec = apply(rna[,5:16], 1, sd)
heatscatter(sd_vec, mean_vec)
abline(v=2, col='red')
abline(h=-4, col='red')
dev.off()


pdf('sd_vs_mean.pdf')
mean_vec = apply(rna[,5:16], 1, mean)
sd_vec = apply(rna[,5:16], 1, sd)
heatscatter(sd_vec, mean_vec)
abline(v=2, col='red')
abline(h=-4, col='red')
dev.off()


for (chr in rownames(table(rna[,1]))){
#chr='chr19'
t = which(rna[,1]==chr)
x = as.matrix(rna[t, 5:16])

tss = rna[t, 2]
m = apply(x,1,mean)
s = apply(x,1,sd)

tt3 = which(m>-4 & s>2)
tt4 = which(m>-4 & s<=2)
tt2 = which(m<=-4 & s>2)
tt1 = which(m<=-4 & s<=2)

used_id_1 = m<=-4 & s<=2
used_id_2 = m<=-4 & s>2
used_id_3 = m>-4 & s>2
used_id_4 = m>-4 & s<=2

used_id_mat = cbind(used_id_1, used_id_2, used_id_3, used_id_4)



	print(paste('chr:', chr))
	for (g in c(1:4)){
		#g=1
		print(paste('group:', g))
		used_id_g = used_id_mat[,g]
		for (c in c(1:12)){
			print(paste('cell:', c))
			dataset_name = paste('vision_rna.', chr, '.', g, '.', c, '.Rdat', sep='')
			output_name = paste('vision_rna.', chr, '.', g, '.', c, '.gene_ccRE.txt', sep='')
			load(dataset_name)
			rna_g = rna[used_id_g,]
			print(summary(rt))
			info_all = cbind(rt$pair, rt$sel)
			extract_gene_ccRE = function(x){
				gene_tss = as.matrix(rna_g[x[3],c(1:4)])
				gene_ccRE = x[2]*200
				info_pos = c(gene_tss, gene_ccRE, x[3], x[4], x[5])
				#print(info_pos)
				return(info_pos)
			}
			info_pos_mat_g = t(apply(info_all, 1, function(x) extract_gene_ccRE(x)))
			write.table(info_pos_mat_g, output_name, quote=FALSE, col.names=FALSE, row.names=FALSE, sep='\t')
		}
	}
}


