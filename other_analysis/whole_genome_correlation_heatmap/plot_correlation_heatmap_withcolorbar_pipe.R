####################################################
### get parameters
args = commandArgs(trailingOnly=TRUE)
signal_file_list = args[1]
signal_high_color = args[2]
signal_low_color = args[3]
heatmap_boarder_col = args[4]
hclust_method = args[5]
correlation_method = args[6]
output_filename = args[7]

#signal_file_list = 'pkn_list.txt'
#signal_high_color = 'red'
#signal_low_color = 'white'
#heatmap_boarder_col = 'black'
#hclust_method = 'average'
#correlation_method = 'pearson'
#output_filename = 'test_cor_heat.png'

#Rscript plot_correlation_heatmap_withcolorbar.R pkn_list.txt red white black average pearson test_cor_heat.png


####################################################
### use rect to plot heatmaps
color_heatmap = function(color_matrix, outputname, format, border_color){
	format(outputname, width = 35, height = 35) ### output name
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

####################################################
### get color bar based on mark
get_color_bar_mark = function(mk){
	if (mk == 'h3k4me3'){
		color_bar_mark = rgb(255/255,0/255,0/255)
	} else if (mk == 'h3k4me1') {
		color_bar_mark = rgb(250/255,250/255,0/255)
	} else if (mk == 'h3k27ac') {
		color_bar_mark = rgb(250/255,150/255,0/255)
	} else if (mk == 'h3k27me3') {
		color_bar_mark = rgb(0/255,0/255,225/255)
	} else if (mk == 'h3k36me3') {
		color_bar_mark = rgb(0/255,150/255,0/255)
	} else if (mk == 'h3k9me3') {
		color_bar_mark = rgb(100/255,100/255,100/255)
	} else if (mk == 'atac') {
		color_bar_mark = rgb(200/255,50/255,150/255)
	} else if (mk == 'dnase') {
		color_bar_mark = rgb(0/255,200/255,200/255)
	} else if (mk == 'ctcf') {
		color_bar_mark = rgb(200/255,0/255,250/255)
	} else if (mk == 'wgbs') {
		color_bar_mark = rgb(30/255,144/255,255/255)
	} else {
		color_bar_mark = rgb(0/255,0/255,0/255)
	}
	return(color_bar_mark)
}

####################################################
### get color bar based on cell type
get_color_bar_ct = function(ct){
	if (ct == 'LSK'){
		color_bar_ct = rgb(34/255,139/255,34/255)
	} else if (ct == 'HPC7') {
		color_bar_ct = rgb(255/255,127/255,0/255)
	} else if (ct == 'CMP') {
		color_bar_ct = rgb(238/255,118/255,0/255)
	} else if (ct == 'MEP') {
		color_bar_ct = rgb(205/255,102/255,0/255)
	} else if (ct == 'G1E') {
		color_bar_ct = rgb(238/255,99/255,99/255)
	} else if (ct == 'ER4') {
		color_bar_ct = rgb(205/255,85/255,85/255)
	} else if (ct == 'CFUE') {
		color_bar_ct = rgb(255/255,48/255,48/255)
	} else if (ct == 'ERY_fl') {
		color_bar_ct = rgb(238/255,44/255,44/255)
	} else if (ct == 'ERY') {
		color_bar_ct = rgb(255/255,0/255,0/255)
	} else if (ct == 'CFUMK') {
		color_bar_ct = rgb(139/255,115/255,85/255)
	} else if (ct == 'iMK') {
		color_bar_ct = rgb(139/255,90/255,43/255)
	} else if (ct == 'MK_fl') {
		color_bar_ct = rgb(139/255,69/255,19/255)
	} else if (ct == 'GMP') {
		color_bar_ct = rgb(0/255,139/255,139/255)
	} else if (ct == 'CLP') {
		color_bar_ct = rgb(147/255,112/255,219/255)
	} else if (ct == 'MON') {
		color_bar_ct = rgb(24/255,116/255,205/255)
	} else if (ct == 'NEU') {
		color_bar_ct = rgb(79/255,148/255,205/255)
	} else if (ct == 'B') {
		color_bar_ct = rgb(139/255,28/255,98/255)
	} else if (ct == 'NK') {
		color_bar_ct = rgb(139/255,0/255,139/255)
	} else if (ct == 'T_CD4') {
		color_bar_ct = rgb(122/255,55/255,139/255)
	} else if (ct == 'T_CD8') {
		color_bar_ct = rgb(104/255,34/255,139/255)
	} 
	return(color_bar_ct)
}

####################################################
### get color bar based on mark
change_ct_names = function(ct){
	if (ct == 'B_SPL') {
		color_bar_mark = 'B'
	} else if (ct == 'CFU_E_ad') {
		color_bar_mark = 'CFUE'
	} else if (ct == 'CFUMK') {
		color_bar_mark = 'CFUMK'
	} else if (ct == 'CLP') {
		color_bar_mark = 'CLP'
	} else if (ct == 'CMP') {
		color_bar_mark = 'CMP'
	} else if (ct == 'ER4') {
		color_bar_mark = 'ER4'
	} else if (ct == 'ERY_ad') {
		color_bar_mark = 'ERY'
	} else if (ct == 'ERY_fl') {
		color_bar_mark = 'ERY_fl'
	} else if (ct == 'G1E') {
		color_bar_mark = 'G1E'
	} else if (ct == 'GMP') {
		color_bar_mark = 'GMP'
	} else if (ct == 'HPC7') {
		color_bar_mark = 'HPC7'
	} else if (ct == 'LSK_BM') {
		color_bar_mark = 'LSK'
	} else if (ct == 'MEP') {
		color_bar_mark = 'MEP'
	} else if (ct == 'MK_imm_ad') {
		color_bar_mark = 'iMK'
	} else if (ct == 'MK_mat_fl') {
		color_bar_mark = 'MK_fl'
	} else if (ct == 'MONO_BM') {
		color_bar_mark = 'MON'
	} else if (ct == 'NEU') {
		color_bar_mark = 'NEU'
	} else if (ct == 'NK_SPL') {
		color_bar_mark = 'NK'
	} else if (ct == 'T_CD4_SPL') {
		color_bar_mark = 'T_CD4'
	} else if (ct == 'T_CD8_SPL') {
		color_bar_mark = 'T_CD8'
	} 
	return(color_bar_mark)
}

###### read file list file
data = read.table(signal_file_list, header = F)
data_matrix = NULL
data_name = c()
color_bar_mark = c()
color_bar_ct = c()
ct_name_vec = c()
for (i in c(1: dim(data)[1])){
        file = paste( toString(data[i,1]), sep='')
        print(file)
        print(i)
        ###### read signal track
        d_tmp = scan(file)
        data_matrix = cbind(data_matrix, d_tmp)
        ###### split file name
        filename_split_vec = unlist(strsplit(file, "[.]"))
        ct_tmp = filename_split_vec[1]
        ct_tmp_modified = change_ct_names(ct_tmp)
        mk_tmp = filename_split_vec[2]
        data_name[i] = paste(ct_tmp_modified, mk_tmp, sep='_')
        ##### get color bar based on mark
        color_bar_mark = rbind(color_bar_mark, get_color_bar_mark(mk_tmp))
        color_bar_ct = rbind(color_bar_ct, get_color_bar_ct(ct_tmp_modified))
        ct_name_vec = rbind(ct_name_vec, paste(ct_tmp_modified, mk_tmp, sep='_'))
}


colnames(data_matrix) = data_name
print(dim(data_matrix))
set.seed(2017)
#used_id = sample(dim(data_matrix)[1], 10000)

###### get correlation matrix between samples
cor_matrix = cor(data_matrix, method = correlation_method)
###### get distance matrix between samples from correlation matrix
dist_cor = as.dist(1 - cor_matrix)
###### cluster samples 
hclust_cor = hclust(dist_cor, method = hclust_method)
new_ct_order = hclust_cor$order


### filter range
filter_range = max(cor_matrix) - min(cor_matrix)
### subtract min(filter matrix) (Entropy smaller -> less white filter)
filter_percent = (cor_matrix - min(cor_matrix) ) / filter_range
filter_percent_label = (cbind(seq(min(cor_matrix), max(cor_matrix), length.out=100), seq(min(cor_matrix), max(cor_matrix), length.out=100)) - min(cor_matrix) ) / filter_range

####################################################
###### convert signal matrix to color matrix (With filter)
###### get r channel score
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

###### reorder correlation heatmap by the correlation-dist  
signal_matrix_color_reorder = signal_matrix_color[new_ct_order, new_ct_order]
color_bar_mark_reorder = color_bar_mark[new_ct_order]
color_bar_ct_reorder = color_bar_ct[new_ct_order]

###### merge corlor bar and correlation matrix
signal_matrix_color_reorder = cbind(signal_matrix_color_reorder, color_bar_mark_reorder, color_bar_ct_reorder)

###### add colnames row names
colnames(signal_matrix_color_reorder) = c(ct_name_vec[new_ct_order], 'mark', 'cell-type')
ct_name_vec_reorder = ct_name_vec[new_ct_order]
rownames(signal_matrix_color_reorder) = ct_name_vec_reorder

rownames(signal_matrix_color_key) = signal_matrix_color_key[,1]

####################################################
###### plot heatmap
heatmap_save_type = pdf
color_heatmap(signal_matrix_color_reorder, output_filename, heatmap_save_type, heatmap_boarder_col)

color_heatmap(signal_matrix_color_key, paste(output_filename,'.colorkey.pdf', sep=''), heatmap_save_type, heatmap_boarder_col)







