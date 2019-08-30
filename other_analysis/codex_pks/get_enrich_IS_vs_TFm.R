library(pheatmap)
library(RColorBrewer)

### get tf chip-seq pk total number
#tf_files = read.table('../file_list.sub.txt', header=F)
#tf_files_names = apply(((tf_files)), 1, function(x) unlist(strsplit(x, '.bb'))[1])
#tf_files_names = apply(((tf_files)), 1, function(x) as.character(x))


JI_mat = c()
union_mat = c()
tf_files_names_exist = c()
file_used = c()
for (i in 0:256){
### input name
input_file_name = paste('JI_matrix_folder_IS_vs_TFmodule/JI_matrix_IS_vs_TFmodule.', i, '.txt', sep='')
### read intersect count of each pk
if(file.exists(input_file_name)){
d = read.table(input_file_name, header=F)
d_JI = (d[,3])
d_union = (d[,2])
JI_mat = rbind(JI_mat, d_JI)
union_mat = rbind(union_mat, d_union)
tf_files_names_exist = c(tf_files_names_exist, i)
file_used = c(file_used, TRUE)
}else{
file_used = c(file_used, FALSE)
}
}


JI_mat = as.matrix(JI_mat)
rownames(JI_mat) = 0:256
colnames(JI_mat) = 1:20

dir.create('IS_vs_TFmodule_JI_plots')
### plot hist JI
pdf('IS_vs_TFmodule_JI_plots/tf.JI_mat.hist.IS_vs_TFmodule.pdf', width=14)
par(mfrow=c(1,2))
hist(JI_mat, breaks=50)
hist(log10(JI_mat+1e-4), breaks=50)
dev.off()


### plot pheatmap
corlo_lim_scale_max = max(JI_mat[-c(256,257),]+1e-4)
corlo_lim_scale_min = min(JI_mat+1e-4)
breaksList_scale = seq(corlo_lim_scale_min, corlo_lim_scale_max, by = 0.001)
my_colorbar_scale=colorRampPalette(c('white', 'blue'))(n = length(breaksList_scale))
JI_dist = dist(1-t(JI_mat))
#JI_dist = dist(t(log10(JI_mat+1e-4)))
pdf('IS_vs_TFmodule_JI_plots/tf.JI_mat.pheatmap.IS_vs_TFmodule.pdf', width = 14, height = 30)
pheatmap((JI_mat+1e-4), clustering_distance_cols=JI_dist, treeheight_col=0, color=my_colorbar_scale, breaks = breaksList_scale, cluster_rows=FALSE, cluster_cols=TRUE, annotation_names_col = FALSE, show_colnames = FALSE)
dev.off()

pdf('IS_vs_TFmodule_JI_plots/tf.JI_mat.pheatmap.IS_vs_TFmodule.withname.pdf', width = 14, height = 30)
pheatmap((JI_mat+1e-4), clustering_distance_cols=JI_dist, treeheight_col=0, color=my_colorbar_scale, breaks = breaksList_scale, cluster_rows=FALSE, cluster_cols=TRUE, annotation_names_col = FALSE, show_colnames = TRUE)
dev.off()

cluster_colors = c('#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', '#46f0f0', '#f032e6', '#bcf60c', '#fabebe', '#008080', '#e6beff', '#9a6324', '#fffac8', '#800000', '#aaffc3', '#808000', '#ffd8b1', '#000075', '#808080')#, '#ffffff', '#000000')
pdf('IS_vs_TFmodule_JI_plots/tf.JI_mat.pheatmap.IS_vs_TFmodule.colcol.pdf', width = 14, height = 5)
pheatmap(as.matrix(rbind(hclust(JI_dist)$order, hclust(JI_dist)$order)), annotation_names_col = FALSE, color = cluster_colors, cluster_rows=FALSE, cluster_cols=FALSE)
dev.off()



corlo_lim_scale_max = max(log10(JI_mat+1e-4))
corlo_lim_scale_min = min(log10(JI_mat+1e-4))
breaksList_scale = seq(corlo_lim_scale_min, corlo_lim_scale_max, by = 0.001)
my_colorbar_scale=colorRampPalette(c('white', 'blue'))(n = length(breaksList_scale))
pdf('IS_vs_TFmodule_JI_plots/tf.JI_mat.pheatmap.log10.IS_vs_TFmodule.pdf', width = 14, height = 30)
pheatmap(log10(JI_mat+1e-4), clustering_distance_cols=JI_dist, color=my_colorbar_scale, breaks = breaksList_scale, cluster_rows=FALSE, cluster_cols=FALSE, annotation_names_col = FALSE)
dev.off()


'''
### log10 scale
corlo_lim_scale_max = max(log10(JI_mat+1e-4))
corlo_lim_scale_min = min(log10(JI_mat+1e-4))
breaksList_scale = seq(corlo_lim_scale_min, corlo_lim_scale_max, by = 0.001)
my_colorbar_scale=colorRampPalette(c('white', 'blue'))(n = length(breaksList_scale))

JI_dist = dist(t(log10(1-JI_mat)+0.001))
pdf('tf_JI_plots/tf.JI_mat.pheatmap.log10.pdf', width = 30, height = 30)
pheatmap(log10(JI_mat+1e-4), clustering_distance_cols=JI_dist, clustering_distance_rows=JI_dist, color=my_colorbar_scale, breaks = breaksList_scale, annotation_names_col = FALSE)
dev.off()
'''



