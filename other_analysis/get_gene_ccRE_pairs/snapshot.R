#snapshot.R
library(pheatmap)
library(MASS)


d = read.table('rnaTPM.txt', header=TRUE, sep=' ')

d0 = d
#d = d[,c(1,2,3,4,11,7,12,5,8,16,15,6,10,9,13,14)]

d = d[,c(1,2,3,4,11,7,12,5,8,6,10,9,13,14)]

pdf('TMP_hist.pdf')
hist(as.matrix(d[,-c(1:4)]), breaks=50)
abline(v=-1, col='red', lwd=1.5, lty=2)
box()
dev.off()

d_binary = (as.matrix(d[,-c(1:4)]) > -1 )*1

d_binary_index = apply(d_binary, 1, function(x) paste(x, collapse='_'))
d_binary_index_0 = d_binary_index
length(table(d_binary_index))

pdf('TMP_index_hist.pdf')
hist(log10(table(d_binary_index)), breaks=50)
abline(v=log10(100), col='red', lwd=1.5, lty=2)
box()
dev.off()

index_set_lim = 100
rare_index = rownames(table(d_binary_index))[table(d_binary_index)<index_set_lim]

for (index in rare_index){
	d_binary_index[d_binary_index_0==index] = 'X_X_X'
}

my_colorbar=colorRampPalette(c('white', 'red'))(n = 128)
pdf(paste('TMP_snapshot.norescuing.pdf', sep=''))
pheatmap(as.matrix(d[order(d_binary_index),-c(1:4)]), color=my_colorbar, cluster_rows = FALSE, cluster_cols = FALSE,annotation_names_row=FALSE,annotation_names_col=FALSE,show_rownames=FALSE,show_colnames=FALSE)
dev.off()


### qda rescuing
d_sig = d[,-c(1:4)]
d_sig = d_sig + matrix(rnorm(dim(d_sig)[1]*dim(d_sig)[2]), dim(d_sig)[1],dim(d_sig)[2])

label_vec_factor = factor(d_binary_index)
qda_input_mat = as.data.frame(d_sig)
qda_input_mat = cbind(qda_input_mat, label_vec_factor)
snapshot_qda_class_mat = c()
for (i in c(1:200)){
#snapshot_qda.fit = qda(label_vec_factor ~ CFU_E_ad+CFUMK+CMP+ERY_fl+GMP+MK_imm_ad+LSK_BM+MEP+MONO_BM+NEU+ER4+G1E, data = qda_input_mat)
snapshot_qda.fit = qda(label_vec_factor ~ CFU_E_ad+CFUMK+CMP+ERY_fl+GMP+MK_imm_ad+LSK_BM+MEP+MONO_BM+NEU, data = qda_input_mat)
snapshot_qda.class = predict(snapshot_qda.fit, qda_input_mat)$class
label_vec_factor = snapshot_qda.class

if (i%%20==0){
snapshot_qda_class_mat = cbind(snapshot_qda_class_mat, table(snapshot_qda.class))
print(table(snapshot_qda.class))
rare_clust = rownames(table(label_vec_factor))[table(label_vec_factor)<index_set_lim]
for (c in rare_clust){
label_vec_factor = as.matrix(label_vec_factor)
label_vec_factor[label_vec_factor==c] = 'X_X_X'
label_vec_factor = factor(label_vec_factor)
}
}

qda_input_mat = as.data.frame(d_sig)
qda_input_mat = cbind(qda_input_mat, label_vec_factor)
}




my_colorbar=colorRampPalette(c('white', 'red'))(n = 128)
pdf(paste('TMP_snapshot.pdf', sep=''))
pheatmap(as.matrix(d[order(label_vec_factor),-c(1:4)]), color=my_colorbar, cluster_rows = FALSE, cluster_cols = FALSE,annotation_names_row=FALSE,annotation_names_col=FALSE,show_rownames=FALSE,show_colnames=FALSE)
dev.off()

output_mat = cbind(d[order(label_vec_factor),c(1:4)], label_vec_factor[order(label_vec_factor)], d[order(label_vec_factor),-c(1:4)])
colnames(output_mat)[5] = 'c'
write.table(output_mat, 'rnaTPM.snapshot.txt', col.names=TRUE, row.names=FALSE, quote=FALSE, sep='\t' )

