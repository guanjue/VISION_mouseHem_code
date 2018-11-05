
### read index set matrix
read_color = function(x){
	rgb_color_int = as.numeric(unlist(strsplit(x, ',')))
	rgb_color = rgb(rgb_color_int[1],rgb_color_int[2],rgb_color_int[3],max=255)
	return(rgb_color)
}


d = read.table('pknorm_2_16lim_ref1mo_0424_lesshet.state', sep=' ')
d = read.table('IDEAS_1_16.wg.state', sep=' ')

d_state_mat = d[,-c(1,2,3,4,25)]

print(sum(rowSums(d_state_mat)>0) / dim(d_state_mat)[1])

colnames(d_state_mat) = c('LSK', 'HPC7', 'CMP', 'MEP', 'G1E', 'ER4', 'CFUE', 'ERY', 'ERY_fl', 'CFUMK', 'iMK', 'MK_fl', 'GMP', 'MON', 'NEU', 'CLP', 'NK', 'B', 'T_CD4', 'T_CD8')

head(d_state_mat)

d_state_mat_non0 = d_state_mat
d_state_mat_non0[d_state_mat!=0]=1

d_state_mat_non0_count = colSums(d_state_mat_non0)
d_state_mat_non0_count_percent = d_state_mat_non0_count / dim(d_state_mat)[1]

ideas_state_color = 'ideas_range_color.txt'
ideas_state_color = 'ideas_range_color.txt'

print('set heatmap colors')
rgb_col_num = read.table(ideas_state_color,header=F)
rgb_col_num = rgb_col_num[,3]
rgb_col_num = rgb_col_num[c((length(rgb_col_num)-1):1,length(rgb_col_num))]
rgb_col_num = as.matrix(rgb_col_num)
#print(rgb_col_num)
rgb_col=apply(rgb_col_num,1,function(x) read_color(x))
my_colorbar=colorRampPalette(rgb_col)(n = dim(rgb_col_num)[1])


ideas_state_matrix_uniq_sort = table(d_state_mat[,1])

counts_matrix = c()
counts_index_matrix = c()
for (i in c(1: dim(d_state_mat)[2]) ){
	print(i)
	### extract ith cell type data
	ideas_state_matrix_table_tmp = as.matrix(d_state_mat[,i])

	table_tmp = c()
	for (j in c( 1: length(ideas_state_matrix_uniq_sort)) ){
		print(j)
		### count the number of cREs have jth IDEAS state
		table_tmp[j] = sum(ideas_state_matrix_table_tmp==(j-1))
	}

	### vector to matrix
	counts_matrix = rbind(counts_matrix, table_tmp)

}

counts_matrix_t_sort = t(counts_matrix)
counts_matrix_t_sort = counts_matrix_t_sort/dim(d_state_mat)[1]
colnames(counts_matrix_t_sort) = colnames(d_state_mat)
rownames(counts_matrix_t_sort) = c( 1: length(ideas_state_matrix_uniq_sort))-1

counts_matrix_t_sort_non0_sum = colSums(counts_matrix_t_sort>0)
counts_matrix_t_sort_0_sum = colSums(counts_matrix_t_sort==0)

pdf('whole_genome_state_count.pdf', dim(d_state_mat)[2]/2+0.5, dim(d_state_mat)[2]/2+0.5)
bb = barplot(counts_matrix_t_sort, col=my_colorbar, ylim=c(0,1), xaxt='n')
text(bb,counts_matrix_t_sort[1,]-0.08,round(colSums(counts_matrix_t_sort[-1,]),3),cex=2, srt = 90)
axis(side=1,at=bb, labels=colnames(counts_matrix_t_sort), las = 2, font=2)
dev.off()






counts_matrix = c()
counts_index_matrix = c()
for (i in c(1: dim(d_state_mat)[2]) ){
	print(i)
	### extract ith cell type data
	ideas_state_matrix_table_tmp = as.matrix(d_state_mat[,i])

	table_tmp = c()
	for (j in c( 2: length(ideas_state_matrix_uniq_sort)) ){
		print(j)
		### count the number of cREs have jth IDEAS state
		table_tmp[j] = sum(ideas_state_matrix_table_tmp==(j-1))
	}

	### vector to matrix
	counts_matrix = rbind(counts_matrix, table_tmp)

}

counts_matrix_t_sort = t(counts_matrix)
counts_matrix_t_sort = counts_matrix_t_sort/dim(d_state_mat)[1]
colnames(counts_matrix_t_sort) = colnames(d_state_mat)
rownames(counts_matrix_t_sort) = c( 1: length(ideas_state_matrix_uniq_sort))-1

counts_matrix_t_sort_non0_sum = colSums(counts_matrix_t_sort>0)
counts_matrix_t_sort_0_sum = colSums(counts_matrix_t_sort==0)



pdf('whole_genome_state_count_notall0.pdf', dim(d_state_mat)[2]/2+0.5, dim(d_state_mat)[2]/2+0.5)
bb = barplot(counts_matrix_t_sort[-1,], col=my_colorbar[-1], ylim=c(0,0.2), xaxt='n')
text(bb,counts_matrix_t_sort[1,]-0.08,round(colSums(counts_matrix_t_sort[-1,]),3), cex=2, srt = 90)
axis(side=1,at=bb, labels=colnames(counts_matrix_t_sort), las=2, font=2)
dev.off()





