### get parameters
args = commandArgs(trailingOnly=TRUE)
index_matrix_ideas_state_inputfile = args[1]
index_matrix_index_inputfile = args[2]
signal_input_list = args[3]
ideas_state_color = args[4]
cREs_IDEASpro_outfile = args[5]
siglim = as.numeric(args[6])
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
class(index_matrix) = 'numeric'
pdf('signal_hist.pdf')
hist(index_matrix, breaks=100)
dev.off()
#siglim = 5
barplot_lim = 100000
index_matrix[index_matrix<=siglim] = 0
index_matrix[index_matrix>siglim] = 1

non_ccRE = sum(rowSums(index_matrix)==0)
print('non_ccREs')
print(non_ccRE)
print('non_ccREs prop')
print(non_ccRE/dim(index_matrix)[1])


print(dim(index_matrix))
print('read ideas state matrix file')
signal_matrix_od = as.matrix(read.table(index_matrix_ideas_state_inputfile, header=FALSE))
signal_matrix_od = signal_matrix_od[,c(-(12+4),-(16+4))]
### extract signal matrix without info
signal_matrix = signal_matrix_od[ , c(5:dim(signal_matrix_od)[2]) ]
### convert to numeric matrix
class(signal_matrix) = 'numeric'

print(head(signal_matrix))
print('repress state')
print(sum((signal_matrix==c(0,2,20,3,16,22))*1))
print(dim(index_matrix))
print(dim(signal_matrix))
#index_matrix[signal_matrix==0] = 0
#index_matrix[signal_matrix==2] = 0
#index_matrix[signal_matrix==20] = 0
#index_matrix[signal_matrix==3] = 0
#index_matrix[signal_matrix==16] = 0
#index_matrix[signal_matrix==22] = 0
#index_matrix[signal_matrix==4] = 0
#index_matrix[signal_matrix==1] = 0

print(dim(signal_matrix))
###### read colnames file
print('signal list')
colname_file = read.table(signal_input_list, header=F)
print(colname_file)
colname = colname_file[,2]
colname = c('LSK','HPC7','CMP','MEP','G1E','ER4','CFUE','ERY','ERY_fl','CFUMK','iMK','GMP','MON','NEU','NK','B','TCD4','TCD8')
colnames(signal_matrix) = colname

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
	ideas_state_matrix_table_tmp = ideas_state_matrix_table_tmp[ideas_state_matrix_table_tmp!=0]
	table_tmp = c()
	for (j in c( 1: length(ideas_state_matrix_uniq_sort)) ){
		### count the number of cREs have jth IDEAS state
		table_tmp[j] = sum(ideas_state_matrix_table_tmp==(j-1))
	}

	### vector to matrix
	counts_matrix = rbind(counts_matrix, table_tmp)

	### extract ith cell type data
	print('atac-pk start')
	ideas_state_matrix_table_tmp = as.matrix(signal_matrix[,i])
	index_state_matrix_table_tmp = as.matrix(index_matrix[,i])
	
	index_state_matrix_table_tmp = index_state_matrix_table_tmp[ideas_state_matrix_table_tmp!=0]
	index_state_matrix_table_tmp = index_state_matrix_table_tmp[index_state_matrix_table_tmp!=0]

	table_tmp = length(index_state_matrix_table_tmp)

	### vector to matrix
	counts_index_matrix = rbind(counts_index_matrix, table_tmp)

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


print(counts_index_matrix_t)
### save figure
pdf(paste(cREs_IDEASpro_outfile,'.index.s3normpk.', toString(siglim), '.pdf',sep=''), width=dim(signal_matrix)[2]*1/20+3, height=dim(signal_matrix)[2]*1/20+3)
barplot(counts_index_matrix_t, col='deepskyblue1', ylim=c(0,barplot_lim),las=2)
#box()
dev.off()


### save figure
pdf(paste(cREs_IDEASpro_outfile, '.ideas.s3normpk.', toString(siglim), '.pdf', sep=''), width=dim(signal_matrix)[2]*1/20+3, height=dim(signal_matrix)[2]*1/20+3)
barplot(counts_matrix_t, col=my_colorbar, ylim=c(0,barplot_lim),las=2)
#box()
dev.off()

#counts_matrix_t_sort = counts_matrix_t[,order(counts_index_matrix_t, decreasing=TRUE)]
#png(paste(cREs_IDEASpro_outfile, '.reorder.s3normpk.', toString(siglim), '.png', sep=''), dim(signal_matrix)[2]*80+5, dim(signal_matrix)[2]*80+5)
#barplot(counts_matrix_t_sort, col=my_colorbar, ylim=c(0,barplot_lim))
#dev.off()


#counts_matrix_t_sort_non0_sum = colSums(counts_matrix_t_sort>0)
#counts_matrix_t_sort_0_sum = colSums(counts_matrix_t_sort==0)
#pdf(paste(cREs_IDEASpro_outfile, '.reorder.s3normpk.', toString(siglim), '.pdf', sep=''), dim(signal_matrix)[2]*1.5+0.5, dim(signal_matrix)[2]*1.5+0.5)
#bb = barplot(counts_matrix_t_sort, col=my_colorbar, ylim=c(0,barplot_lim))
#text(bb,counts_matrix_t_sort[1,]-2000,colSums(counts_matrix_t_sort[-1,]),cex=1.5)
#dev.off()

### transpose matrix
#counts_index_matrix_t_sort = as.matrix(counts_index_matrix_t[order(counts_index_matrix_t, decreasing=TRUE)])
#names = colnames(signal_matrix)
#names_sort = names[order(counts_index_matrix_t, decreasing=TRUE)]
#rownames(counts_index_matrix_t_sort) = names_sort
#counts_index_matrix_t_sort = t(counts_index_matrix_t_sort)

#print(counts_index_matrix_t_sort)
### save figure
#png(paste(cREs_IDEASpro_outfile, '.reorder.index.s3normpk.', toString(siglim), '.png', sep=''), dim(signal_matrix)[2]*80+5, dim(signal_matrix)[2]*80+5)
#barplot(counts_index_matrix_t_sort, col='deepskyblue1', ylim=c(0,barplot_lim))
#dev.off()

#pdf(paste(cREs_IDEASpro_outfile, '.reorder.index.s3normpk.', toString(siglim), '.pdf', sep=''), dim(signal_matrix)[2]*1.5+0.5, dim(signal_matrix)[2]*1.5+0.5)
#bb = barplot(counts_index_matrix_t_sort, col='deepskyblue1', ylim=c(0,barplot_lim))
#text(bb,counts_index_matrix_t_sort+2000,counts_index_matrix_t_sort,cex=1.5)
#dev.off()

# Rscript plot_ct_IDEASpro_bar.R atac_20cell.bed.ideas.matrix.txt.freq.indexed.sort.txt ideas_list.txt ideas_range_color.txt bar.pdf

