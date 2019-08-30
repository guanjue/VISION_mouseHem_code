library(pheatmap)
library(RColorBrewer)

### get tf chip-seq pk total number
tf_files = read.table('../file_list.txt', header=F)
#tf_files_names = apply(((tf_files)), 1, function(x) unlist(strsplit(x, '.bb'))[1])
tf_files_names = apply(((tf_files)), 1, function(x) as.character(x))


JI_mat = c()
union_mat = c()
tf_files_names_exist = c()
file_used = c()
for (i in tf_files_names){
i = as.character(i)
print(i)
### input name
input_file_name = paste('JI_matrix_folder/JI_matrix_', i, '.txt', sep='')
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

JI_mat = JI_mat[,file_used]
JI_mat = as.matrix(JI_mat)
rownames(JI_mat) = tf_files_names_exist
colnames(JI_mat) = tf_files_names_exist

dir.create('tf_JI_plots')
### plot hist JI
pdf('tf_JI_plots/tf.JI_mat.hist.800.pdf', width=14)
par(mfrow=c(1,2))
hist(JI_mat, breaks=50)
hist(log10(JI_mat+1e-4), breaks=50)
dev.off()


### get JI distance
JI_dist = dist(1-(JI_mat))
hc = hclust(JI_dist, "complete")
###
dist_within = c()
nclustr = 50
for (i in c(1:nclustr)){
cluster_id = cutree(hc, i)
dist_within_i = 0
for (j in c(1:i)){
dist_within_j = sum(1-JI_mat[cluster_id==j, cluster_id==j])
dist_within_i = dist_within_i + dist_within_j
}
dist_within[i] = dist_within_i
}
### plot JI dist curve
pdf('tf_JI_plots/tf.JI_dist.curve.800.pdf')
par(mfrow=c(1,1))
plot(1:nclustr, dist_within, log='')
lines(1:nclustr, dist_within)
dev.off()

############
### cut tree
ncluster = 30
cluster_id = cutree(hc, ncluster)

### plot pheatmap
corlo_lim_scale_max = max(JI_mat+1e-4)
corlo_lim_scale_min = min(JI_mat+1e-4)
breaksList_scale = seq(corlo_lim_scale_min, corlo_lim_scale_max, by = 0.001)
my_colorbar_scale=colorRampPalette(c('white', 'blue'))(n = length(breaksList_scale))
JI_dist = dist(1-(JI_mat))
#JI_dist = dist(t(log10(JI_mat+1e-4)))
pdf('tf_JI_plots/tf.JI_mat.pheatmap.800.pdf', width = 100, height = 100)
pheatmap((JI_mat+1e-4), clustering_distance_cols=JI_dist, clustering_distance_rows=JI_dist, color=my_colorbar_scale, breaks = breaksList_scale, annotation_names_col = FALSE)
dev.off()

tf_files_names_exist_len = apply(cbind(tf_files_names_exist), 1, function(x) length(unlist(strsplit(x,''))))
tf_files_names_max = tf_files_names_exist[tf_files_names_exist_len==max(tf_files_names_exist_len)]
label_matrix = cbind(cluster_id, cluster_id)
colnames(label_matrix) = c(tf_files_names_max,tf_files_names_max)
pdf('tf_JI_plots/tf.JI_mat.pheatmap.cluster.800.pdf', width = 7, height = 100)
pheatmap(label_matrix, clustering_distance_rows=JI_dist, annotation_names_col = FALSE, color = colorRampPalette(rev(brewer.pal(n = ncluster, name ="Paired")))(ncluster))
dev.off()






### log10 scale
corlo_lim_scale_max = max(log10(JI_mat+1e-4))
corlo_lim_scale_min = min(log10(JI_mat+1e-4))
breaksList_scale = seq(corlo_lim_scale_min, corlo_lim_scale_max, by = 0.001)
my_colorbar_scale=colorRampPalette(c('white', 'blue'))(n = length(breaksList_scale))

JI_dist = dist(t(log10(1-JI_mat)+1e-4))
JI_dist = dist(1-(JI_mat))

pdf('tf_JI_plots/tf.JI_mat.pheatmap.log10.800.pdf', width = 30, height = 30)
pheatmap(log10(JI_mat+1e-4), clustering_distance_cols=JI_dist, clustering_distance_rows=JI_dist, color=my_colorbar_scale, breaks = breaksList_scale, annotation_names_col = FALSE)
dev.off()




