d = read.table('pknorm_2_16lim_ref1mo_0424_lesshet.para')

mark_colnames = c('atac', 'ctcf', 'h3k27ac', 'h3k27me3', 'h3k36me3', 'h3k4me1', 'h3k4me3', 'h3k9me3')
sig_matrix = as.matrix(d[,c(2:9)] / d[,1])

library(pheatmap)
my_colorbar=colorRampPalette(c('white', 'blue'))(n = 128)
col_breaks = c(seq(0, 2000,length=33))

rownames(sig_matrix) = c(0:(dim(sig_matrix)[1]-1))
colnames(sig_matrix) = mark_colnames


pdf('ideas_heatmap.pdf')
pheatmap(sig_matrix, color=my_colorbar, cluster_cols = FALSE, cluster_rows=TRUE, annotation_names_row=TRUE, annotation_names_col=TRUE,show_rownames=TRUE,show_colnames=TRUE)
dev.off()


write.table(sig_matrix, file='ideas_state_mean_signal.txt', quote=FALSE, sep='\t', col.names=TRUE, row.names=TRUE)