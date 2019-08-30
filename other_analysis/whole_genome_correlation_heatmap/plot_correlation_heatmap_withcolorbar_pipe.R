####################################################
### get parameters
args = commandArgs(trailingOnly=TRUE)
signal_file_list = args[1]
ct_color_list = args[2]
signal_high_color = args[3]
signal_low_color = args[4]
heatmap_boarder_col = args[5]
hclust_method = args[6]
correlation_method = args[7]
output_filename = args[8]
figure_size = as.numeric(args[9])
sample_size = as.numeric(args[10])

#signal_file_list = 'pkn_list.txt'
#ct_color_list = 'ct_color_list.txt'
#signal_high_color = 'red'
#signal_low_color = 'white'
#heatmap_boarder_col = 'black'
#hclust_method = 'average'
#correlation_method = 'pearson'
#output_filename = 'test_cor_heat.png'
#figure_size = 35
#sample_size = 100000

#Rscript plot_correlation_heatmap_withcolorbar.R file_list.txt ct_color_list.txt red white black average pearson test_0319.pdf 35 100000


####################################################
### use rect to plot heatmaps
color_heatmap = function(color_matrix, outputname, format, border_color, figure_size){
	format(outputname, width = figure_size, height = figure_size) ### output name
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
	mk = tolower(mk)
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


###### read ct color
ct_color_list = read.table(ct_color_list, header = F, sep='\t')
####################################################
### get color bar based on cell type
get_color_bar_ct = function(ct, ct_color_list){
	if (sum(tolower(ct_color_list[,1])==tolower(ct))!=0){
		color0 = ct_color_list[tolower(ct_color_list[,1])==tolower(ct),2]
		color1 = as.numeric( unlist( strsplit(as.character(color0), ',') ) )
		color_bar_ct = rgb(color1[1]/255,color1[2]/255,color1[3]/255)
	} else {
		color_bar_ct = rgb(0/255,0/255,0/255)
	}
	return(color_bar_ct)
}


###### read file list file
data = read.table(signal_file_list, header = F)
data_matrix = NULL
data_name = c()
color_bar_mark = c()
color_bar_ct = c()
ct_name_vec = c()
for (i in c(1: dim(data)[1])){
        file = as.character(data[i,3])
        print(file)
        print(i)
        ###### read signal track
        d_tmp = scan(file)
        data_matrix = cbind(data_matrix, d_tmp)
        ###### get ct name and mk name
        ct_tmp = data[i,1]
        mk_tmp = data[i,2]
        data_name[i] = paste(ct_tmp, mk_tmp, sep='_')
        ##### get color bar based on mark
        color_bar_mark = rbind(color_bar_mark, get_color_bar_mark(mk_tmp))
        color_bar_ct = rbind(color_bar_ct, get_color_bar_ct(ct_tmp, ct_color_list))
        ct_name_vec = rbind(ct_name_vec, paste(ct_tmp, mk_tmp, sep='_'))
}


colnames(data_matrix) = data_name
print(dim(data_matrix))
set.seed(2017)
used_id = sample(dim(data_matrix)[1], sample_size)

###### get correlation matrix between samples
cor_matrix = cor(data_matrix[used_id,], method = correlation_method)
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
color_heatmap(signal_matrix_color_reorder, output_filename, heatmap_save_type, heatmap_boarder_col, figure_size)

color_heatmap(signal_matrix_color_key, paste(output_filename,'.colorkey.pdf', sep=''), heatmap_save_type, heatmap_boarder_col, figure_size)







