
library(pheatmap)

d1 = read.table('~/group/projects/vision/withMETH/run_IDEAS_7ct_nometh_result/run_IDEAS_7ct_nometh.state', header=F, sep=' ')
d2 = read.table('~/group/projects/vision/withMETH/run_IDEAS_7ct_result/run_IDEAS_7ct.state', header=F, sep=' ')

set.seed(2018)
sample_id = sample(dim(d1)[1], 100000)
state1 = d1[sample_id,-c(1:4,dim(d1)[2])]
state2 = d2[sample_id,-c(1:4,dim(d1)[2])]

state_pair_count = function(v1, v2){
	v1_s = 0:32
	v2_s = 0:32
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
	count_mat = matrix(nrow=33, ncol=33)
	n_c = length(v1)
	for (i in 0:32){
		print(i)
		for (j in 0:32){
			print(j)
			#print((sum(v1==i & v2==j)+10))
			#print((sum(v1==i) / n_c * sum(v2==j)+10))
			#print((sum(v1==i & v2==j)+10)/(sum(v1==i) / n_c * sum(v2==j)+10))
			count_mat[as.numeric(i)+1, as.numeric(j)+1] = (sum(v1==i & v2==j)+10) / (sum(v1==i) / n_c * sum(v2==j)+10)
		}
	}
	return(count_mat)
}

LSK_c_mat = state_pair_count(state1[,1], state2[,2])

rownames(LSK_c_mat) = c(0:32)
colnames(LSK_c_mat) = c(0:32)

write.table(LSK_c_mat, 'all_c_mat.txt', sep='\t', quote=F)

all_c_mat = state_pair_count(as.vector(as.matrix(state1[,])), as.vector(as.matrix(state2[,])))

rownames(all_c_mat) = c(0:32)
colnames(all_c_mat) = c(0:32)

write.table(all_c_mat, 'all_c_mat.txt', sep='\t', quote=F)

my_colorbar=colorRampPalette(c('white', 'red'))(n = 128)

LSK_c_mat0 = LSK_c_mat
lim_c =100000
LSK_c_mat0[LSK_c_mat0>lim_c]=lim_c
LSK_c_mat0[is.na(LSK_c_mat0)]=1
#LSK_c_mat0[LSK_c_mat0<1]=1

pdf('LSK_c_mat.pdf', width=9)
pheatmap((LSK_c_mat0), color=my_colorbar, cluster_cols = TRUE,cluster_rows=TRUE,annotation_names_row=TRUE,annotation_names_col=TRUE,show_rownames=TRUE,show_colnames=TRUE,clustering_method = 'complete')
dev.off()


all_c_mat0 = all_c_mat
lim_c =100000
all_c_mat0[all_c_mat0>lim_c]=lim_c
all_c_mat0[is.na(all_c_mat0)]=min(all_c_mat0[!is.na(all_c_mat0)])
#all_c_mat0[all_c_mat0<1]=1


color_n = 50
my_colorbar=colorRampPalette(c('white', 'red'))(n = color_n)

breaksList = seq(0, color_n, by = 1)
pdf('all_c_mat.pdf', width=9)
pheatmap((all_c_mat0), color=my_colorbar, breaks = breaksList, clustering_method = 'complete', cluster_cols = TRUE,cluster_rows=TRUE,annotation_names_row=TRUE,annotation_names_col=TRUE,show_rownames=TRUE,show_colnames=TRUE)
dev.off()


breaksList = seq(-8, 8, by = 0.1)

color_n = length(breaksList)
my_colorbar=colorRampPalette(c('blue', 'white', 'red'))(n = color_n)

pdf('all_c_mat.log2.pdf', width=9)
pheatmap(log2(all_c_mat0), color=my_colorbar, breaks = breaksList, clustering_method = 'complete', cluster_cols = T,cluster_rows=T,annotation_names_row=TRUE,annotation_names_col=TRUE,show_rownames=TRUE,show_colnames=TRUE)
dev.off()

pdf('all_c_mat.log2.nocluster.pdf', width=9)
pheatmap(log2(all_c_mat0), color=my_colorbar, breaks = breaksList, clustering_method = 'complete', cluster_cols = F,cluster_rows=F,annotation_names_row=TRUE,annotation_names_col=TRUE,show_rownames=TRUE,show_colnames=TRUE)
dev.off()


