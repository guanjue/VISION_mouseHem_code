### get parameters
args = commandArgs(trailingOnly=TRUE)

input_matrix = args[1]
output_heatmap = args[2]

data = read.table(input_matrix, header=F, sep='\t')
data_sig = data[,-1]
rownames(data_sig) = data[,1]
colnames(data_sig) = c('CH12', 'Fetal_Liver', 'MEL')

library(pheatmap)

my_colorbar=colorRampPalette(c('white', 'red'))(n = 128)
col_breaks = c(seq(0, 2000,length=33))

png(output_heatmap)
pheatmap((data_sig), color=my_colorbar, cluster_cols = FALSE,cluster_rows=FALSE,annotation_names_row=FALSE,annotation_names_col=TRUE,show_rownames=TRUE,show_colnames=TRUE)
dev.off()

png(paste(output_heatmap, '.cluster.png',sep=''))
pheatmap(log2(data_sig), color=my_colorbar, cluster_cols = FALSE,cluster_rows=TRUE,annotation_names_row=FALSE,annotation_names_col=TRUE,show_rownames=TRUE,show_colnames=TRUE, clustering_method='average')
dev.off()
