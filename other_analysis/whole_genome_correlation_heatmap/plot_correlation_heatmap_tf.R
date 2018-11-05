####################################################
### use rect to plot heatmaps
color_heatmap = function(color_matrix, outputname, format, border_color){
	#format(outputname, width = 35, height = 35) ### output name
	format(outputname, width = 100, height = 100)
	par(mar=c(10,0.5,0.5,10)) ### set heatmap margins
	colbin_len = 10 ### column bin size
	rowbin_len = 10 ### row bin size
	### row reverse
	color_matrix = color_matrix[nrow(color_matrix):1,]
	### plot areas
	plot(c(0, dim(color_matrix)[2]*colbin_len), c(0, dim(color_matrix)[1]*rowbin_len), xaxt = "n", yaxt = "n", xaxs="i", yaxs="i", type = "n", xlab = "", ylab = "",main = "")
	### add color matrix colname as heatmap colname
	axis(1, c(1 : dim(color_matrix)[2])*colbin_len-0.5*colbin_len, colnames(color_matrix), las = 2, col.axis = "black", tick=FALSE)
	axis(4, c(1 : dim(color_matrix)[1])*colbin_len-0.5*colbin_len, rownames(color_matrix), las = 2, col.axis = "black", tick=FALSE)
	### use for loop to add rectangle with different color
	for (coln in c(1 : dim(color_matrix)[2])){ ### loop columns
		for (rown in c(1 : dim(color_matrix)[1])){ ### loop rows
			### add rectangle
			rect( (coln-1)*colbin_len, (rown-1)*rowbin_len, coln*colbin_len, rown*rowbin_len, col = color_matrix[rown, coln], border=border_color, lwd = 0 )
		}
	}
	dev.off()
}

###### read file list file
data = read.table('allconcate.sort.mergebed.mat.txt', header = F)

set.seed(2018)
used_id = sample(dim(data)[1], 100000)

data_r = data[used_id,]
ct = read.table('bb_list_ct.txt', header=F)
tf = read.table('bb_list_tf.txt', header=F)
filename = read.table('bb_list_filename.txt', header=F)

tf_ct_vec = c()
for (i in c(1:dim(ct)[1])){
	tf_ct_vec[i] = paste(tf[i,1], ct[i,1], filename[i,1], sep=':')
}

library(RColorBrewer)
color = grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)]
col_vector_tf = apply(t(col2rgb(sample(color, length(table(tf))))), 1, function(x) rgb(x[1]/255, x[2]/255, x[3]/255))
col_vector_ct = apply(t(col2rgb(sample(color, length(table(ct))))), 1, function(x) rgb(x[1]/255, x[2]/255, x[3]/255))

###get tf label color
tf_table = table(tf)

col_list_tf = list()
for (i in c(1:length(tf_table))){
	key = rownames(tf_table)[i]
	value = col_vector_tf[i]
	col_list_tf[[key]] = value
}

col_vec_tf = c()
for (i in c(1:dim(tf)[1])){
	tf_i = tf[i,1]
	col_vec_tf[i] = unlist(col_list_tf[tf_i])
}

###get ct label color
ct_table = table(ct)

col_list_ct = list()
for (i in c(1:length(ct_table))){
	key = rownames(ct_table)[i]
	value = col_vector_ct[i]
	col_list_ct[[key]] = value
}

col_vec_ct = c()
for (i in c(1:dim(ct)[1])){
	ct_i = ct[i,1]
	col_vec_ct[i] = unlist(col_list_ct[ct_i])
}



data_matrix=as.matrix(data_r[,-c(1:3)])
correlation_method = 'pearson'

cor_matrix = cor(data_matrix, method = correlation_method)
###### get distance matrix between samples from correlation matrix
dist_cor = as.dist(1 - cor_matrix)
###### cluster samples 
hclust_method = 'average'
hclust_cor = hclust(dist_cor, method = hclust_method)
new_ct_order = hclust_cor$order


### filter range
filter_range = max(cor_matrix) - min(cor_matrix)
### subtract min(filter matrix) (Entropy smaller -> less white filter)
filter_percent = (cor_matrix - min(cor_matrix) ) / filter_range
filter_percent_label = (cbind(seq(min(cor_matrix), max(cor_matrix), length.out=100), seq(min(cor_matrix), max(cor_matrix), length.out=100)) - min(cor_matrix) ) / filter_range

signal_high_color = 'red'
signal_low_color = 'white'

signal_high_color_rgb = colorRamp(signal_high_color)(1)
signal_low_color_rgb = colorRamp(signal_low_color)(1)

### signal to color
rh = signal_high_color_rgb[1]/255
rl = signal_low_color_rgb[1]/255
### add filter: od_color * filter + bg_color * (1 - filter)
r = rh * filter_percent + rl * (1-filter_percent)
r_label = rh * filter_percent_label + rl * (1-filter_percent_label)
### set upper & lower limit for r
r[r>1] = 1
r[r<0] = 0
###### get g channel score
### signal to color
gh = signal_high_color_rgb[2]/255
gl = signal_low_color_rgb[2]/255
### add filter: od_color * (1-filter) + bg_color * filter
g = gh * filter_percent + gl * (1-filter_percent)
g_label = gh * filter_percent_label + gl * (1-filter_percent_label)
### set upper & lower limit for g
g[g>1] = 1
g[g<0] = 0
###### get b channel score
### signal to color
bh = signal_high_color_rgb[3]/255
bl = signal_low_color_rgb[3]/255
### add filter: od_color * (1-filter) + bg_color * filter
b = bh * filter_percent + bl * (1-filter_percent)
b_label = bh * filter_percent_label + bl * (1-filter_percent_label)
### set upper & lower limit for b
b[b>1] = 1
b[b<0] = 0
###### creat rgb color matrix
### for heatmap
signal_matrix_color = NULL
for (i in seq(1,dim(r)[2])){
	#print(i)
	### convert each r,g,b vector to a rgb matrix column
	cor_col_tmp = rgb( r[,i], g[,i], b[,i] )
	### cbind column to a rgb matrix
	signal_matrix_color = cbind(signal_matrix_color, cor_col_tmp)
}
### for color key
signal_matrix_color_key = NULL
for (i in seq(1,dim(r_label)[2])){
	### convert each r,g,b vector to a rgb matrix column
	cor_col_tmp = rgb( r_label[,i], g_label[,i], b_label[,i] )
	### cbind column to a rgb matrix
	signal_matrix_color_key = cbind(signal_matrix_color_key, cor_col_tmp)
}



signal_matrix_color_reorder = signal_matrix_color[new_ct_order,new_ct_order]
color_bar_mark_reorder = col_vec_tf[new_ct_order]
color_bar_ct_reorder = col_vec_ct[new_ct_order]

###### merge corlor bar and correlation matrix
signal_matrix_color_reorder = cbind(signal_matrix_color_reorder, color_bar_mark_reorder, color_bar_ct_reorder)
colnames(signal_matrix_color_reorder) = c(tf_ct_vec[new_ct_order], 'mark', 'cell-type')
rownames(signal_matrix_color_reorder) = tf_ct_vec[new_ct_order]


heatmap_save_type = pdf
heatmap_boarder_col = NA
color_heatmap(signal_matrix_color_reorder, 'codex_tf_cor_tfmerged.colorbar.pdf', heatmap_save_type, heatmap_boarder_col)

signal_matrix_color_reorder_noct = cbind(signal_matrix_color[new_ct_order,new_ct_order], color_bar_mark_reorder)
colnames(signal_matrix_color_reorder_noct) = c(tf_ct_vec[new_ct_order], 'mark')
rownames(signal_matrix_color_reorder_noct) = tf_ct_vec[new_ct_order]

color_heatmap(signal_matrix_color_reorder_noct, 'codex_tf_cor_tfmerged.noct.colorbar.pdf', heatmap_save_type, heatmap_boarder_col)

signal_matrix_color_reorder_03 = signal_matrix_color[new_ct_order,new_ct_order]
signal_matrix_color_reorder_03[cor_matrix[new_ct_order,new_ct_order]<0.3] = "#FFFFFF"
signal_matrix_color_reorder_03 = cbind(signal_matrix_color_reorder_03, color_bar_mark_reorder, color_bar_ct_reorder)
colnames(signal_matrix_color_reorder_03) = c(tf_ct_vec[new_ct_order], 'mark', 'cell-type')
rownames(signal_matrix_color_reorder_03) = tf_ct_vec[new_ct_order]

color_heatmap(signal_matrix_color_reorder_03, 'codex_tf_cor_tfmerged03.colorbar.pdf', heatmap_save_type, heatmap_boarder_col)


signal_matrix_color_reorder_noct_03 = signal_matrix_color[new_ct_order,new_ct_order]
signal_matrix_color_reorder_noct_03[cor_matrix[new_ct_order,new_ct_order]<0.3] = "#FFFFFF"
signal_matrix_color_reorder_noct_03 = cbind(signal_matrix_color_reorder_noct_03, color_bar_mark_reorder)
colnames(signal_matrix_color_reorder_noct_03) = c(tf_ct_vec[new_ct_order], 'mark')
rownames(signal_matrix_color_reorder_noct_03) = tf_ct_vec[new_ct_order]
color_heatmap(signal_matrix_color_reorder_noct_03, 'codex_tf_cor_tfmerged03.noct.colorbar.pdf', heatmap_save_type, heatmap_boarder_col)














###### cluster samples
library(pheatmap)
colnames(data_matrix) = tf_ct_vec
colnames(cor_matrix) = tf_ct_vec
rownames(cor_matrix) = tf_ct_vec

my_colorbar=colorRampPalette(c('white', 'red'))(n = 128)
col_breaks = c(seq(0, 2000,length=33))

pdf('codex_tf_cor.pdf', width=100, height=100)
pheatmap(cor_matrix, clustering_distance_rows=dist_cor, clustering_distance_cols=dist_cor, color=my_colorbar, annotation_names_col = FALSE)
dev.off()

png('codex_tf_cor.png', width=7000, height=7000)
pheatmap(cor_matrix, clustering_distance_rows=dist_cor, clustering_distance_cols=dist_cor, color=my_colorbar, annotation_names_col = FALSE)
dev.off()

pdf('cor_hist.pdf')
hist(cor_matrix, breaks=50)
dev.off()


png('codex_tf_cor_03.png', width=7000, height=7000)
cor_matrix_03 = cor_matrix
cor_matrix_03[cor_matrix<0.3]=0
pheatmap(cor_matrix_03, clustering_distance_rows=dist_cor, clustering_distance_cols=dist_cor, color=my_colorbar, annotation_names_col = FALSE)
dev.off()

pdf('codex_tf_cor_03.pdf', width=100, height=100)
pheatmap(cor_matrix_03, clustering_distance_rows=dist_cor, clustering_distance_cols=dist_cor, color=my_colorbar, annotation_names_col = FALSE)
dev.off()


tf_name = colnames(table(tf))

data_matrix_ctmerge = c()
for (i in c(1:length(tf_name))){
	print(tf_name[i])
	if (sum(tf == tf_name[i]) > 1){
		data_matrix_tmp = (rowSums(data_matrix[,tf==tf_name[i]]) > 0)*1
	} else{
		data_matrix_tmp = data_matrix[,tf==tf_name[i]]
	}
	data_matrix_ctmerge = cbind(data_matrix_ctmerge, data_matrix_tmp)
}

colnames(data_matrix_ctmerge) = tf_name

cor_matrix_ctmerge = cor(data_matrix_ctmerge, method = correlation_method)
###### get distance matrix between samples from correlation matrix
dist_cor_ctmerge = as.dist(1 - cor_matrix_ctmerge)





heatmap_save_type = pdf
heatmap_boarder_col = NA
color_heatmap(signal_matrix_color_reorder, 'codex_tf_cor_tfmerged.colorbar.pdf', heatmap_save_type, heatmap_boarder_col)



color_heatmap(signal_matrix_color_reorder, 'codex_tf_cor_tfmerged.noct.colorbar.pdf', heatmap_save_type, heatmap_boarder_col)


color_heatmap(signal_matrix_color_key, paste(output_filename,'.colorkey.pdf', sep=''), heatmap_save_type, heatmap_boarder_col)



































































pdf('codex_tf_cor_tfmerged.pdf', width=100, height=100)
pheatmap(cor_matrix_ctmerge, clustering_distance_rows=dist_cor_ctmerge, clustering_distance_cols=dist_cor_ctmerge, color=my_colorbar, annotation_names_col = FALSE)
dev.off()


png('codex_tf_cor_tfmerged.png', width=2500, height=2500)
pheatmap(cor_matrix_ctmerge, clustering_distance_rows=dist_cor_ctmerge, clustering_distance_cols=dist_cor_ctmerge, color=my_colorbar, annotation_names_col = FALSE)
dev.off()

pdf('cor_hist_tfmerge.pdf')
hist(cor_matrix, breaks=50)
dev.off()


png('codex_tf_cor_03_tfmerged.png', width=2500, height=2500)
cor_matrix_03 = cor_matrix
cor_matrix_ctmerge[cor_matrix_ctmerge<0.3]=0
pheatmap(cor_matrix_ctmerge, clustering_distance_rows=dist_cor_ctmerge, clustering_distance_cols=dist_cor_ctmerge, color=my_colorbar, annotation_names_col = FALSE)
dev.off()

pdf('codex_tf_cor_03_tfmerged.pdf', width=30, height=30)
pheatmap(cor_matrix_ctmerge, clustering_distance_rows=dist_cor_ctmerge, clustering_distance_cols=dist_cor_ctmerge, color=my_colorbar, annotation_names_col = FALSE)
dev.off()


