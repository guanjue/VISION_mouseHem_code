eRP0 = read.table('eRP0.txt', header=F)[,2]

tpm = read.table('gene_tss.2MB.IDEAS_states.tpm.txt', header=F)
ideas_state = read.table('gene_tss.IDEAS_states.2MB.IDEAS_state.txt', header=F)

get_eRP_mat = function(state_vec){
	t = apply(as.matrix(state_vec),2,function(x) eRP0[x+1])
	return(t)
}

n = dim(ideas_state)[1]

ideas_state_eRPs = t(apply(ideas_state[1:n,], 1, get_eRP_mat))

tpm_eRPs = cbind(tpm[1:n,], ideas_state_eRPs)

tpm_eRPs_cor = apply(tpm_eRPs, 1, function(x) cor(x[1:12], x[13:24]))

write.table(tpm_eRPs_cor, 'tpm_eRPs_cor.txt', quote=F, sep='\t', col.names=F, row.names=F)

pdf('hist_cor.0.pdf')
hist(tpm_eRPs_cor, breaks=30)
abline(v=0.2, col='red')
box()
dev.off()

pdf('hist_cor.abs.0.pdf')
hist(abs(tpm_eRPs_cor), breaks=50)
abline(v=0.2, col='red')
box()
dev.off()


gene_ccRE_table = read.table('gene_ccRE_pair_table.0.txt', header=F)
ccRE_pre_gene = length(table(gene_ccRE_table[,2]))/length(table(gene_ccRE_table[,1]))
print(ccRE_pre_gene)

c = table(gene_ccRE_table[,1])
summary(cbind(c))


