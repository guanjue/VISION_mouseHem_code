library(rGREAT)
library(qvalue)

bed_list = as.matrix(read.table('bed_list.txt', header = F))
print(bed_list)

great_MGIPhenotype = c()
great_MGIPhenotype_hyper_fc = c()
great_MGIPhenotype_binom_fc = c()
great_MGIPhenotype_hyper_p = c()

for (i in seq(1, dim(bed_list)[1])){
	bed_file = bed_list[i]
	print(bed_file)
	### read bed file
	bed = read.table(bed_file, header = F)[,c(1,2,3)]
	### submit to GREAT
	job = submitGreatJob(bed, species = 'mm10', request_interval = 300, rule='basalPlusExt', adv_upstream=5.0, adv_downstream=1.0, adv_span=100)
	#pdf(paste(toString(bed_file), '.dist.pdf', sep=''), width=15, height=5)
	#par(mfrow = c(1, 3))
	#res = plotRegionGeneAssociationGraphs(job)
	#dev.off()
	### get GO analysis result
	tb = getEnrichmentTables(job, ontology = 'Mouse Phenotype') 
	### GO Molecular Function
	tb_x = tb[['Mouse Phenotype']]
	tb_x_padj = as.matrix(tb_x[order(tb_x[,1]),][8])
	#tb_x_padj = qvalue(p = as.matrix(tb_x[order(tb_x[,1]),][8]))$qvalues
	tb_x_hyper_fc = as.matrix(tb_x[order(tb_x[,1]),][12])
	tb_x_binom_fc = as.matrix(tb_x[order(tb_x[,1]),][6])
	tb_x_hyper_p = as.matrix(tb_x[order(tb_x[,1]),][15])
	### cbind 
	great_MGIPhenotype = cbind(great_MGIPhenotype, tb_x_padj)
	great_MGIPhenotype_hyper_fc = cbind(great_MGIPhenotype_hyper_fc, tb_x_hyper_fc)
	great_MGIPhenotype_hyper_p = cbind(great_MGIPhenotype_hyper_p, tb_x_hyper_p)
	great_MGIPhenotype_binom_fc = cbind(great_MGIPhenotype_binom_fc, tb_x_binom_fc)
}

tb_x = tb[['Mouse Phenotype']]
id = tb_x[order(tb_x[,1]),1]
names = tb_x[order(tb_x[,1]),2]
total_num = tb_x[order(tb_x[,1]),9]
great_MGIPhenotype = cbind(id, names, total_num, great_MGIPhenotype)
great_MGIPhenotype_hyper_fc = cbind(id, names, total_num, great_MGIPhenotype_hyper_fc)
great_MGIPhenotype_hyper_p = cbind(id, names, total_num, great_MGIPhenotype_hyper_p)
great_MGIPhenotype_binom_fc = cbind(id, names, total_num, great_MGIPhenotype_binom_fc)

write.table(great_MGIPhenotype, file = 'Mouse_Phenotype.txt', sep='\t', quote=F, col.names=F, row.names=F)

write.table(great_MGIPhenotype_hyper_fc, file = 'great_MGIPhenotype_hyper_fc.txt', sep='\t', quote=F, col.names=F, row.names=F)

write.table(great_MGIPhenotype_hyper_p, file = 'great_MGIPhenotype_hyper_p.txt', sep='\t', quote=F, col.names=F, row.names=F)

wrwrite.table(great_MGIPhenotype_binom_fc, file = 'great_MGIPhenotype_binom_fc.txt', sep='\t', quote=F, col.names=F, row.names=F)

ite.table(total_num, file = 'Mouse_Phenotype_totalnum.txt', sep='\t', quote=F, col.names=F, row.names=F)


data=read.table('great_MGIPhenotype_binom_fc.txt', header = F, sep='\t', quote='')
bp_p = (as.matrix(data[, c(-1,-2,-3)]))
go_total_num = data[,3]
rownames(bp_p) = data[,2]
colnames(bp_p) = seq(0,dim(bp_p)[2]-1)
#bp_p = bp_p[,-dim(bp_p)[2]]
rowmax = apply(bp_p,1,max)
bp_p_sig = bp_p[as.logical((rowmax >= 0) * (go_total_num<=100000)),]

bp_p_sig_lim = bp_p_sig
bp_p_sig_lim[bp_p_sig_lim>100] =100

library(pheatmap)
my_colorbar = colorRampPalette(c('white', 'red'))(n = 128)
png('binom_fc.png', width=9000,height=9000)
pheatmap(bp_p_sig_lim, color=my_colorbar, cluster_cols = F,cluster_rows=T,annotation_names_row=FALSE,annotation_names_col=F,show_rownames=T,show_colnames=T, clustering_method='ward.D2')
dev.off()

pdf('binom_fc.pdf', width=50,height=1000)
pheatmap(bp_p_sig_lim, color=my_colorbar, cluster_cols = F,cluster_rows=T,annotation_names_row=FALSE,annotation_names_col=F,show_rownames=T,show_colnames=T, clustering_method='ward.D2')
dev.off()



data=read.table('Mouse_Phenotype.txt', header = F, sep='\t', quote='')
bp_p_fdr = apply(data[, c(-1,-2,-3)], 2, function(x) p.adjust(x, method = 'fdr'))
#bp_p = -log10(as.matrix(data[, c(-1,-2,-3)]))
bp_p = -log10(as.matrix(bp_p_fdr))
go_total_num = data[,3]
rownames(bp_p) = data[,2]
colnames(bp_p) = seq(0,dim(bp_p)[2]-1)
#bp_p = bp_p[,-dim(bp_p)[2]]
rowmax = apply(bp_p,1,max)
bp_p_sig = bp_p[as.logical((rowmax >= 0) * (go_total_num<=100000)),]

bp_p_sig_lim = bp_p_sig
bp_p_sig_lim[bp_p_sig_lim>100] =100

library(pheatmap)
my_colorbar = colorRampPalette(c('white', 'red'))(n = 128)
png('binom_p.png', width=9000,height=9000)
pheatmap(bp_p_sig_lim, color=my_colorbar, cluster_cols = F,cluster_rows=T,annotation_names_row=FALSE,annotation_names_col=F,show_rownames=T,show_colnames=T, clustering_method='ward.D2')
dev.off()

pdf('binom_p.pdf', width=50,height=1000)
pheatmap(bp_p_sig_lim, color=my_colorbar, cluster_cols = F,cluster_rows=T,annotation_names_row=FALSE,annotation_names_col=F,show_rownames=T,show_colnames=T, clustering_method='ward.D2')
dev.off()


bp_p_sig_lim  = bp_p_sig_lim
rowmax_bp = apply(bp_p_sig_lim, 1, max)
bp_p_sig_lim_min001 = bp_p_sig_lim[(rowmax_bp)>=0,]
colnames(bp_p_sig_lim_min001) = bed_list
rownames(bp_p_fdr) = data[,2]
write.table(bp_p_sig_lim_min001, file = 'Mouse_Phenotype_max_all.txt', sep='\t', quote=F, col.names=T, row.names=T)

bp_p_sig_lim  = bp_p_sig_lim
rowmax_bp = apply(bp_p_sig_lim, 1, max)
bp_p_sig_lim_min001 = bp_p_sig_lim[(rowmax_bp)>=2,]
colnames(bp_p_sig_lim_min001) = bed_list
rownames(bp_p_fdr) = data[,2]
write.table(bp_p_sig_lim_min001, file = 'Mouse_Phenotype_max_001.txt', sep='\t', quote=F, col.names=T, row.names=T)



data=read.table('great_MGIPhenotype_hyper_p.txt', header = F, sep='\t', quote='')
bp_p = -log10(as.matrix(data[, c(-1,-2,-3)]))
go_total_num = data[,3]
rownames(bp_p) = data[,2]
colnames(bp_p) = seq(0,dim(bp_p)[2]-1)
#bp_p = bp_p[,-dim(bp_p)[2]]
rowmax = apply(bp_p,1,max)
bp_p_sig = bp_p[as.logical((rowmax >= 0) * (go_total_num<=100000)),]

bp_p_sig_lim = bp_p_sig
bp_p_sig_lim[bp_p_sig_lim>100] =100

library(pheatmap)
my_colorbar = colorRampPalette(c('white', 'red'))(n = 128)
png('hyper_p.png', width=9000,height=9000)
pheatmap(bp_p_sig_lim, color=my_colorbar, cluster_cols = F,cluster_rows=T,annotation_names_row=FALSE,annotation_names_col=F,show_rownames=T,show_colnames=T, clustering_method='ward.D2')
dev.off()

pdf('hyper_p.pdf', width=50,height=1000)
pheatmap(bp_p_sig_lim, color=my_colorbar, cluster_cols = F,cluster_rows=T,annotation_names_row=FALSE,annotation_names_col=F,show_rownames=T,show_colnames=T, clustering_method='ward.D2')
dev.off()


data=read.table('great_MGIPhenotype_hyper_fc.txt', header = F, sep='\t', quote='')
bp_p = (as.matrix(data[, c(-1,-2,-3)]))
go_total_num = data[,3]
rownames(bp_p) = data[,2]
colnames(bp_p) = seq(0,dim(bp_p)[2]-1)
#bp_p = bp_p[,-dim(bp_p)[2]]
rowmax = apply(bp_p,1,max)
bp_p_sig = bp_p[as.logical((rowmax >= 0) * (go_total_num<=100000)),]

bp_p_sig_lim = bp_p_sig
bp_p_sig_lim[bp_p_sig_lim>100] =100

library(pheatmap)
my_colorbar = colorRampPalette(c('white', 'red'))(n = 128)
png('hyper_fc.png', width=9000,height=9000)
pheatmap(bp_p_sig_lim, color=my_colorbar, cluster_cols = F,cluster_rows=T,annotation_names_row=FALSE,annotation_names_col=F,show_rownames=T,show_colnames=T, clustering_method='ward.D2')
dev.off()

pdf('hyper_fc.pdf', width=50,height=1000)
pheatmap(bp_p_sig_lim, color=my_colorbar, cluster_cols = F,cluster_rows=T,annotation_names_row=FALSE,annotation_names_col=F,show_rownames=T,show_colnames=T, clustering_method='ward.D2')
dev.off()




