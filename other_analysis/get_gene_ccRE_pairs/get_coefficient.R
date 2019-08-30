setwd('/gpfs/group/yzz2/default/legacy/group/projects/vision/rna/get_ccRE_gene_pair_nofilter/')
tpm = read.table('gene_tss.IDEAS_states.tpm.txt', header=F)
ideas_state = read.table('gene_tss.IDEAS_states.IDEAS_state.txt', header=F)

tpm_vec = matrix(as.matrix(tpm), nrow=dim(tpm)[1]*dim(tpm)[2],byrow=T)
ideas_state_vec = matrix(as.matrix(ideas_state), nrow=dim(ideas_state)[1]*dim(ideas_state)[2],byrow=T)

state_type = as.numeric(rownames(table(ideas_state_vec)))

state2matrix = function(x){
	seq0 = rep(0, 27)
	seq0[x+1] = 1
	return(seq0)
}

ideas_state_mat = t(apply(ideas_state_vec, 1, state2matrix))
y_ideas_state = cbind(tpm_vec, ideas_state_mat)
colnames(y_ideas_state) = c('tpm', c(0:26))


y_ideas_state = as.data.frame(y_ideas_state)

y_ideas_state_cor = cor(y_ideas_state)

write.table(cbind(y_ideas_state_cor[-1,1]), 'cor0.txt', quote=F, sep='\t', col.names=F, row.names=T)

y_ideas_state_all = y_ideas_state

y_ideas_state = y_ideas_state[,-2]


fit_lm = lm(tpm ~ .-1, data=y_ideas_state)
summary(fit_lm)

eRP0 = cbind(fit_lm$coefficients)
eRP0 = rbind(0, eRP0)
rownames(eRP0) = c(0:26)

colnames(eRP0) = c('eRP0')

write.table(eRP0, 'eRP0.txt', quote=F, sep='\t', col.names=F, row.names=T)


fit_lm_all = lm(tpm ~ .-1, data=y_ideas_state_all)
summary(fit_lm_all)

eRP_all = cbind(fit_lm_all$coefficients)
rownames(eRP_all) = c(0:26)

colnames(eRP_all) = c('eRP_all')

write.table(eRP_all, 'eRP_all.txt', quote=F, sep='\t', col.names=F, row.names=T)

