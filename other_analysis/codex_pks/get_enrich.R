library(pheatmap)

### get whole genome length
wg = read.table('~/group/genome/mm9/mm9.1to19_X.genome')
wg_len = sum(wg[,2])

### get tf chip-seq pk total number
bg_info = read.table('tf_pk_num_all.txt', header=F)
bg_count = bg_info[,1]
dataset = apply(cbind(as.character(bg_info[,2])), 1, function(x) unlist(strsplit(x, '/'))[2])
dataset1 = apply(cbind(dataset), 1, function(x) unlist(strsplit(x, '.bb'))[1])


enrichment_mat = c()
for (i in c(1:9)){
print(i)
### input name
input_file_name = paste('count_pk_matrix_', i, '.txt', sep='')
### read intersect count of each pk
d = read.table(input_file_name, header=F)
d_len = sum(d[,3]-d[,2])

### read intersect count of pk cluster
fg_count = colSums(d[,-c(1:4)])

### get expect intersect counts
exp_count = d_len / wg_len * bg_count

### get enrichment
enrichment = (fg_count+100) / (exp_count+100)
enrichment_mat = cbind(enrichment_mat, enrichment)
}

enrichment_mat = as.matrix(enrichment_mat)
colnames(enrichment_mat) = c(1:9)
rownames(enrichment_mat) = dataset1

### get enrichment hist
pdf('enrichment_hist.pdf')
hist(rowMeans(log2(enrichment_mat)), breaks=50)
dev.off()

my_colorbar=colorRampPalette(c('white', 'blue'))(n = 128)
col_breaks = c(seq(0, 2000,length=33))

pdf('tf.enrichment.pheatmap.pdf', height = 100)
pheatmap(log2(enrichment_mat), color=my_colorbar, annotation_names_col = FALSE, cluster_cols=FALSE)
dev.off()

pdf('tf.enrichment.pheatmap.lim1.pdf', height = 100)
pheatmap(log2(enrichment_mat[rowMeans(log2(enrichment_mat))>=1,]), color=my_colorbar, annotation_names_col = FALSE, cluster_cols=FALSE)
dev.off()


rmctcf = apply(cbind(dataset1), 1, function(x) !grepl('ctcf', tolower(x)))
enrichment_mat_scale = ((enrichment_mat))
enrichment_mat_scale = t(apply((enrichment_mat), 1, function(x) x-mean(x)))

enrichment_mat_noctcf = enrichment_mat_scale[rmctcf,]

pdf('tf.enrichment.pheatmap.lim1.noctcf.pdf', height = 100)
plot_enrichment_mat_noctcf = (enrichment_mat_noctcf[rowMeans(enrichment_mat_noctcf[,])>0,])
plot_enrichment_mat_noctcf[plot_enrichment_mat_noctcf<0] = 0
pheatmap(plot_enrichment_mat_noctcf, color=my_colorbar, annotation_names_col = FALSE, cluster_cols=FALSE)
dev.off()




