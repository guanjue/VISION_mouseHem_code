data = read.table('pkn_list.txt', header = F)

data_matrix = NULL
data_name = c()
for (i in c(1: dim(data)[1])){
	file = paste( toString(data[i,1]), sep='')
	print(file)
	print(i)
	d_tmp = scan(file)
	data_matrix = cbind(data_matrix, d_tmp)
	data_name[i] = paste(unlist(strsplit(file, "[.]"))[1], unlist(strsplit(file, "[.]"))[2], sep='_')
}

colnames(data_matrix) = data_name

print(dim(data_matrix))

set.seed(2017)
used_id = sample(dim(data_matrix)[1], 500000)

library(pheatmap)
my_colorbar=colorRampPalette(c('white', 'red'))(n = 128)
col_breaks = c(seq(0, 2000,length=33))

cor_matrix = cor(data_matrix[used_id,], method = 'pearson')
pdf('pkn.pc.pdf', width=50, height=50)
pheatmap(cor_matrix, color=my_colorbar, cluster_cols = TRUE,cluster_rows=TRUE,annotation_names_row=TRUE,annotation_names_col=TRUE,show_rownames=TRUE,show_colnames=TRUE)
dev.off()

write.table(as.matrix(cor_matrix), 'pk_cor.matrix.txt', sep='\t', quote = FALSE, row.names = FALSE, col.names = TRUE)

cor_matrix = cor(data_matrix[used_id,], method = 'spearman')
pdf('pkn.sp.pdf', width=50, height=50)
pheatmap(cor_matrix, color=my_colorbar, cluster_cols = TRUE,cluster_rows=TRUE,annotation_names_row=TRUE,annotation_names_col=TRUE,show_rownames=TRUE,show_colnames=TRUE)
dev.off()

write.table(as.matrix(cor_matrix), 'pk_cor.matrix.sp.txt', sep='\t', quote = FALSE, row.names = FALSE, col.names = TRUE)

