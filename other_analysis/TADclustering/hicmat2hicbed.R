d_rc = read.table('/storage/home/g/gzx103/group/HiC/mouse/processed/10kb/G1E-ER4all3.merged.chr19.10kb.matrix', header=F)

bedpe_mat = matrix(data = 0, nrow = (dim(d_rc)[1]*dim(d_rc)[1])/2+dim(d_rc)[1]/2, ncol = 8)

k = 0
for (i in 1:dim(d_rc)[1]){
if (i%%100==0){
print(i)
}
for (j in 1:dim(d_rc)[1]){
if (i >= j){
k = k+1
rc_tmp = d_rc[i,j]
if (rc_tmp>0){
bedpe = c('chr19', i*10000, (i+1)*10000, 'chr19', j*10000, (j+1)*10000, paste(i,j,'_'), rc_tmp)
bedpe_mat[k,] = bedpe
}
}
}
}

write.table(bedpe_mat, 'G1E-ER4all3.merged.chr19.10kb.matrix.bedpe', quote=F, col.names=F, row.names=F, sep='\t')






