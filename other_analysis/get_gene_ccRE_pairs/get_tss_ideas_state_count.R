col_sum_all = rep(0, 27)

for (i in c(1:19)){
print(i)
for (j in c(1:4)){
rdat_name = paste('vision_rna.chr', toString(i), '.', toString(j), '.', '1', '.Rdat', sep='')
load(rdat_name)
col_sum_tmp = colSums(rt$x0)
col_sum_all = col_sum_all + col_sum_tmp
}
print(col_sum_all)
}

write.table(cbind(col_sum_all), 'tss_ideas_state_prop.txt', quote=F, col.names=F, row.names=F, sep='\t')


col_sum_all = rep(0, 27)

for (i in c(1:19)){
print(i)
for (j in c(1:4)){
rdat_name = paste('vision_rna_tss5k.chr', toString(i), '.', toString(j), '.', '1', '.Rdat', sep='')
load(rdat_name)
col_sum_tmp = colSums(rt$x0)
col_sum_all = col_sum_all + col_sum_tmp
}
print(col_sum_all)
}

write.table(cbind(col_sum_all), 'tss5k_ideas_state_prop.txt', quote=F, col.names=F, row.names=F, sep='\t')


col_sum_all = rep(0, 27)

for (i in c(1:19)){
print(i)
for (j in c(1:4)){
rdat_name = paste('vision_rna_tss5k.chr', toString(i), '.', toString(j), '.', '1', '.Rdat', sep='')
load(rdat_name)
col_sum_tmp = colSums(rt$x)
col_sum_all = col_sum_all + col_sum_tmp
}
print(col_sum_all)
}

write.table(cbind(col_sum_all), 'dist5k_ideas_state_prop.txt', quote=F, col.names=F, row.names=F, sep='\t')

rdat_name = paste('vision_rna_tss5k.chr', toString(i), '.', toString(j), '.', '1', '.Rdat', sep='')
load(rdat_name)
d = matrix(NA, nrow = 27, ncol = 27)
for (i in c(1:27)){
	for (j in c(1:27)){
	print(cor(rt$x0[,i],rt$x[,j]))
	d[i,j] = cor(rt$x0[,i],rt$x[,j])
}
}
summary(d)
max(d)
