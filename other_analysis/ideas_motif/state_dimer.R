library(pheatmap)


d0 = read.table('ideas_m2.state', header=T, comment='~')

set.seed(2019)
d = d0[1:1000000,]

dl_allct = d[,-c(1:4, dim(d)[2])]
ct_num = dim(dl_allct)[2]
state_num = length(rownames(table(d[,5])))
ct_names = colnames(d0)[-c(1:4, dim(d)[2])]

state_dimer_allct_mat = matrix(data = 0, nrow = state_num, ncol = state_num)

for (ct_i in 1:ct_num){
print(ct_i)
d_ct = dl_allct[,ct_i]
state_1 = d_ct[-length(d_ct)]
state_2 = d_ct[-1]
for (i in 0:(state_num-1)){
#print(i)
for (j in 0:(state_num-1)){
	state_dimer = sum((state_1==i) & (state_2==j))
	state_dimer_allct_mat[(i+1),(j+1)] = state_dimer_allct_mat[(i+1),(j+1)]+state_dimer
}
}
}


colnames(state_dimer_allct_mat) = c(0:(state_num-1))
rownames(state_dimer_allct_mat) = c(0:(state_num-1))

state_dimer_allct_mat = state_dimer_allct_mat/ct_num


for (ct_i in 1:ct_num){
print(ct_names[ct_i])
dl = dl_allct[,ct_i]
state_1 = dl[-length(dl)]
state_2 = dl[-1]
state_dimer_mat = matrix(data = NA, nrow = state_num, ncol = state_num)
for (i in 0:(state_num-1)){
print(i)
for (j in 0:(state_num-1)){
	state_dimer = sum((state_1==i) & (state_2==j))
	state_dimer_mat[(i+1),(j+1)] = state_dimer
}
}
colnames(state_dimer_mat) = c(0:(state_num-1))
rownames(state_dimer_mat) = c(0:(state_num-1))
my_colorbar=colorRampPalette(c('white', 'blue'))(n = 50)
pdf(paste(ct_names[ct_i], 'test.dimer.pdf', sep=''))
pheatmap(log2(state_dimer_mat+1),cluster_rows=F,cluster_cols=F,color=my_colorbar,show_rownames=TRUE,show_colnames=TRUE,annotation_names_row=FALSE,annotation_names_col=TRUE)
dev.off()
state_dimer_mat_enrich = (state_dimer_mat+1)/(state_dimer_allct_mat+1)
state_dimer_mat_enrich[state_dimer_mat_enrich>10]=10
state_dimer_mat_enrich[1,1]=10
pdf(paste(ct_names[ct_i], 'test.dimer.enrich.pdf', sep=''))
pheatmap((state_dimer_mat_enrich),cluster_rows=F,cluster_cols=F,color=my_colorbar,show_rownames=TRUE,show_colnames=TRUE,annotation_names_row=FALSE,annotation_names_col=TRUE)
dev.off()
}


