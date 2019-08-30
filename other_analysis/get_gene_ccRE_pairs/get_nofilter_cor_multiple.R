eRP0 = read.table('eRP0.txt', header=F)[,2]

get_eRP_mat = function(state_vec){
	t = apply(as.matrix(state_vec),2,function(x) eRP0[x+1])
	return(t)
}


pdf(paste('hist_cor.', 'all', '.pdf', sep=''), height=50)
par(mfrow=c(10,1))
for ( i in c(1:10, 1000)){
print(i)
tpm = read.table(paste('gene_tss.2MB.IDEAS_states.tpm.rand.', i, '.txt', sep=''), header=F)
ideas_state = read.table(paste('gene_tss.IDEAS_states.2MB.IDEAS_state.rand.', i, '.txt', sep=''), header=F)
#
n = dim(ideas_state)[1]
#
ideas_state_eRPs = t(apply(ideas_state[1:n,], 1, get_eRP_mat))
#
tpm_eRPs = cbind(tpm[1:n,], ideas_state_eRPs)
#
tpm_eRPs_cor = apply(tpm_eRPs, 1, function(x) cor(x[1:12], x[13:24]))
write.table(tpm_eRPs_cor, paste('cor.', i, '.txt', sep=''), quote=F, col.names=F, row.names=F, sep='\t')
#
hist(tpm_eRPs_cor, breaks=30, xlim=c(-1,1))
abline(v=0.2, col='red')
box()
}
#
dev.off()


tpm_eRPs_cor1 = read.table('cor.1.txt', header=F)
tpm_eRPs_cor10 = read.table('cor.1000.txt', header=F)

pdf('hist_1kb_bin.pdf')
plot(density(as.matrix(tpm_eRPs_cor1[!is.na(tpm_eRPs_cor1)])), xlim=c(-1,1), ylim=c(0, 1.3), col='blue')
lines(density(as.matrix(tpm_eRPs_cor10[!is.na(tpm_eRPs_cor10)])))
abline(v=0.2, col='red')
dev.off()


pdf('hist_1kb_bin.cdf.pdf')
plot(ecdf(as.matrix(tpm_eRPs_cor1[!is.na(tpm_eRPs_cor1)])), xlim=c(-1,1), ylim=c(0, 1), col='blue')
lines(ecdf(as.matrix(tpm_eRPs_cor10[!is.na(tpm_eRPs_cor10)])))
abline(v=0.2, col='red')
dev.off()
