library(pheatmap)

### get parameters
args = commandArgs(trailingOnly=TRUE)

input_fimo_info = args[1]#'G1Ehigh.G1E.gene.withbinID.tranmat1d.500.txt'
output_heatmap = args[2]#'G1Ehigh.ER4.gene.withbinID.tranmat1d.500.txt'
bin_length = as.numeric(args[3])

d = read.table(input_fimo_info, header=T)
d_max_pos = as.integer((d[,3]+d[,4])/2)

bin_num = dim(d)[1]

fimo_score_mat = matrix(0,bin_num,bin_length/10)

for (i in 1:bin_num){
if (i%%100==0){
print(i)
}
fimo_score_mat[i,as.integer(d_max_pos[i]/10)] = 1
}




my_colorbar=colorRampPalette(c('white', 'black'))(n = 128)
col_breaks = c(seq(0, 2000,length=33))
plot_id_row = seq(1, dim(fimo_score_mat)[1], 1)
plot_id_col = seq(1, dim(fimo_score_mat)[2], 1)
png(output_heatmap)
pheatmap(fimo_score_mat[plot_id_row, plot_id_col], color=my_colorbar, cluster_cols = FALSE,cluster_rows=FALSE, clustering_method = 'average',annotation_names_row=FALSE,annotation_names_col=FALSE,show_rownames=FALSE,show_colnames=FALSE)
dev.off()
