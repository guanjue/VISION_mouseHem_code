d = read.table('ccRE_eRP.table.txt', header=F)
d_sig = d[,-1]
d_sig_r = round(d_sig, 2)
d = cbind(d[,1], d_sig_r)
write.table(d, 'ccRE_eRP.table.r.txt', quote=F, col.names=F, row.names=F, sep='\t')

