
library(pheatmap)

d1 = read.table('~/group/projects/vision/pknorm_2_16lim_ref1mo_sample1mo_0424/pknorm_2_16lim_ref1mo_0424_lesshet.state', header=F, sep=' ')
d2 = read.table('~/group/projects/vision/withMETH/run_IDEAS_result_lesshet/run_IDEAS_lesshet.state', header=F, sep=' ')

set.seed(2018)
sample_id = sample(dim(d1)[1], 100000)
state1 = d1[sample_id,-c(1:4,dim(d1)[2])]
state2 = d2[sample_id,-c(1:4,dim(d1)[2])]

state_pair_count = function(v1, v2){
	v1_s = 0:26
	v2_s = 0:31
	count_mat = matrix(nrow=length(v1_s), ncol=length(v2_s))
	for (i in 1:length(v1_s)){
		print(i)
		for (j in 1:length(v2_s)){
			print(j)
			count_mat[i, j] = sum(v1==i & v2==j)
		}
	}
	return(count_mat)
}

state_pair_count = function(v1, v2){
	count_mat = matrix(nrow=27, ncol=32)
	n_c = length(v1)
	for (i in 0:26){
		print(i)
		for (j in 0:31){
			print(j)
			#print((sum(v1==i & v2==j)+10))
			#print((sum(v1==i) / n_c * sum(v2==j)+10))
			#print((sum(v1==i & v2==j)+10)/(sum(v1==i) / n_c * sum(v2==j)+10))
			count_mat[as.numeric(i)+1, as.numeric(j)+1] = (sum(v1==i & v2==j)+10) / (sum(v1==i) / n_c * sum(v2==j)+10)
		}
	}
	return(count_mat)
}


all_c_mat = state_pair_count(as.vector(as.matrix(state1[,])), as.vector(as.matrix(state2[,])))

rownames(all_c_mat) = c(0:26)
colnames(all_c_mat) = c(0:31)

write.table(all_c_mat, 'all_c_mat_20ct.txt', sep='\t', quote=F)

my_colorbar=colorRampPalette(c('white', 'red'))(n = 128)



all_c_mat0 = all_c_mat
lim_c =100000
all_c_mat0[all_c_mat0>lim_c]=lim_c
all_c_mat0[is.na(all_c_mat0)]=min(all_c_mat0[!is.na(all_c_mat0)])
#all_c_mat0[all_c_mat0<1]=1


color_n = 100
my_colorbar=colorRampPalette(c('white', 'red'))(n = color_n)

breaksList = seq(0, color_n, by = 1)
pdf('all_c_mat_20ct.pdf', width=9)
pheatmap((all_c_mat0), color=my_colorbar, breaks = breaksList, clustering_method = 'ward.D2', cluster_cols = TRUE,cluster_rows=TRUE,annotation_names_row=TRUE,annotation_names_col=TRUE,show_rownames=TRUE,show_colnames=TRUE)
dev.off()




