library(pheatmap)

### get parameters
args = commandArgs(trailingOnly=TRUE)

input_fimo_info = args[1]#'G1Ehigh.G1E.gene.withbinID.tranmat1d.500.txt'
output_heatmap = args[2]#'G1Ehigh.ER4.gene.withbinID.tranmat1d.500.txt'
input_fimo_info = 'cmp9_eryad3.5kbupdown.gata.fimo.info.txt'
output_heatmap = 'cmp9_eryad3.5kbupdown.gata.fimo.heatmap.png'

d = read.table(input_fimo_info, header=T)
bin_length = sum(d[,1]==d[1,1])/2

bin_num = length(unique(d[,1]))


d_fimo_score = -log10(d[,4])
d_fimo_score = (d[,4])

fimo_score_mat = matrix(0,bin_num,bin_length)

for (i in 1:bin_num){
if (i%%100==0){
print(i)
}
n_start = (i-1)*bin_length*2+1
n_end = (i)*bin_length*2
used_id_pos = seq(n_start,n_end,2)
#used_id_neg = seq(n_start,n_end,2)+1
fimo_score_pos = d_fimo_score[used_id_pos]
max_pos = which.max(fimo_score_pos)
scores = rep(0, bin_length)
scores[max_pos] = 10
#fimo_score_neg = d_fimo_score[used_id_neg]
#fimo_score_max = scores#apply(cbind(fimo_score_pos, fimo_score_neg), 1, max)
fimo_score_mat[i,] = scores
}

my_colorbar=colorRampPalette(c('white', 'black'))(n = 128)
col_breaks = c(seq(0, 2000,length=33))
plot_id_row = seq(1, dim(fimo_score_mat)[1], 1)
plot_id_col = seq(1, dim(fimo_score_mat)[2], 1)
png(output_heatmap)
pheatmap(fimo_score_mat[plot_id_row, plot_id_col], color=my_colorbar, cluster_cols = FALSE,cluster_rows=FALSE, clustering_method = 'average',annotation_names_row=FALSE,annotation_names_col=FALSE,show_rownames=FALSE,show_colnames=FALSE)
dev.off()


png('test.png')
plot(1:993, cbind(colMeans(fimo_score_mat)))
dev.off()


