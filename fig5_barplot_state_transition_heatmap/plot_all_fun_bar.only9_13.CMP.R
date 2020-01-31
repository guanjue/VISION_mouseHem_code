### get parameters
args = commandArgs(trailingOnly=TRUE)
index_matrix_ideas_state_inputfile = args[1]
index_matrix_index_inputfile = args[2]
signal_input_list = args[3]
ideas_state_color = args[4]
cREs_IDEASpro_outfile = args[5]

index_matrix_ideas_state_inputfile = 'atac_pk_no0.atac.ideas.txt'
index_matrix_index_inputfile = 'atac_pk_no0.atac.index.txt'
signal_input_list = 'signal_list_atac.txt'
ideas_state_color = 'ideas_range_color.txt'
cREs_IDEASpro_outfile = 'atac_18cell.allfunbar.CMP'

### read index set matrix
read_color = function(x){
	rgb_color_int = as.numeric(unlist(strsplit(x, ',')))
	rgb_color = rgb(rgb_color_int[1],rgb_color_int[2],rgb_color_int[3],max=255)
	return(rgb_color)
}


#################################################### 
############ read input files
####################################################
### read signal matrix file
print('read signal matrix file')
index_matrix_od = as.matrix(read.table(index_matrix_index_inputfile, header=FALSE))
index_matrix = index_matrix_od[ , c(5:dim(index_matrix_od)[2]) ]
print(dim(index_matrix))
print('read ideas state matrix file')
signal_matrix_od = as.matrix(read.table(index_matrix_ideas_state_inputfile, header=FALSE))
signal_matrix_od = signal_matrix_od[,c(-(12+4),-(16+4))]
### extract signal matrix without info
signal_matrix = signal_matrix_od[ , c(5:dim(signal_matrix_od)[2]) ]
print(dim(signal_matrix))
### convert to numeric matrix
class(signal_matrix) = 'numeric'
###### read colnames file
print('signal list')
colname_file = read.table(signal_input_list, header=F)
print(colname_file)
colname = colname_file[,2]
colname = c('LSK','HPC7','CMP','MEP','G1E','ER4','CFUE','ERY','ERY_fl','CFUMK','iMK','GMP','MON','NEU','NK','B','TCD4','TCD8')
colnames(signal_matrix) = colname


used_id = apply(signal_matrix, 1, function(x) 9 %in% x)
index_matrix = index_matrix[used_id,]
signal_matrix = signal_matrix[used_id,]


### index set
print('get uniq index set id')
index_set_id = signal_matrix_od[,4]
index_set_id_uniq = unique(index_set_id)
### sort index
index_set_id_uniq_sort = sort(index_set_id_uniq)

### get unique elements from the matrix
print('get uniq ideas label id')
ideas_state_matrix_flatten = as.vector(signal_matrix)
ideas_state_matrix_uniq = unique(ideas_state_matrix_flatten)
### sort index
ideas_state_matrix_uniq_sort = sort(ideas_state_matrix_uniq)

### set heatmap colors
print('set heatmap colors')
rgb_col_num = read.table(ideas_state_color,header=F)
rgb_col_num = rgb_col_num[,3]
rgb_col_num = rgb_col_num[c((length(rgb_col_num)-1):1,length(rgb_col_num))]
rgb_col_num = as.matrix(rgb_col_num)
#print(rgb_col_num)
rgb_col=apply(rgb_col_num,1,function(x) read_color(x))
my_colorbar=colorRampPalette(rgb_col)(n = dim(rgb_col_num)[1])


###### extract counts matrix
counts_matrix = c()
counts_index_matrix = c()
for (i in c(1: dim(signal_matrix)[2]) ){
### extract ith cell type data
ideas_state_matrix_table_tmp = as.matrix(signal_matrix[,i])
index_state_matrix_table_tmp = as.matrix(index_matrix[,i])

ideas_state_matrix_table_tmp = ideas_state_matrix_table_tmp[index_state_matrix_table_tmp!=0]
#ideas_state_matrix_table_tmp = ideas_state_matrix_table_tmp[ideas_state_matrix_table_tmp!=0]
table_tmp = c()
for (j in c( 1: length(ideas_state_matrix_uniq_sort)) ){
#for (j in c(9+1,13+1) ){
### count the number of cREs have jth IDEAS state
#table_tmp[j] = sum(ideas_state_matrix_table_tmp==(j-1))
table_tmp = c(table_tmp, sum(ideas_state_matrix_table_tmp==(j-1)))
}
### vector to matrix
counts_matrix = rbind(counts_matrix, table_tmp)
### extract ith cell type data
print('atac-pk start')
ideas_state_matrix_table_tmp = as.matrix(signal_matrix[,i])
index_state_matrix_table_tmp = as.matrix(index_matrix[,i])
index_state_matrix_table_tmp = index_state_matrix_table_tmp#[ideas_state_matrix_table_tmp!=0]
index_state_matrix_table_tmp = index_state_matrix_table_tmp#[index_state_matrix_table_tmp!=0]
#table_tmp = length(index_state_matrix_table_tmp)
### vector to matrix
counts_index_matrix = rbind(counts_index_matrix, table_tmp[10])
}

### transpose matrix
counts_matrix_t = t( counts_matrix )
### add colnames
colnames(counts_matrix_t) = colnames(signal_matrix)

### save figure
png(cREs_IDEASpro_outfile, dim(signal_matrix)[2]*100+5, dim(signal_matrix)[2]*100+5)
barplot(counts_matrix_t, col=my_colorbar)
dev.off()

### transpose matrix
counts_index_matrix_t = t( counts_index_matrix )
### add colnames
colnames(counts_index_matrix_t) = colnames(signal_matrix)

barplot_lim = 60000
print(counts_index_matrix_t)
### save figure

counts_matrix_t_sort = counts_matrix_t[,order(counts_index_matrix_t, decreasing=TRUE)]
print(counts_matrix_t_sort)
counts_matrix_t_sort_cumsum = apply(counts_matrix_t_sort, 2, cumsum)
pdf(paste(cREs_IDEASpro_outfile, '.reorder.9_13.pdf', sep=''), height=3.5,width=5)
bb = barplot(counts_matrix_t_sort, col=my_colorbar, ylim=c(0,barplot_lim),las=2)
box()
#text(bb,apply(counts_matrix_t_sort, 2, cumsum),counts_matrix_t_sort,cex=1.5)
dev.off()

# Rscript plot_ct_IDEASpro_bar.R atac_20cell.bed.ideas.matrix.txt.freq.indexed.sort.txt ideas_list.txt ideas_range_color.txt bar.pdf

