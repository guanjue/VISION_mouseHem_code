setwd('/gpfs/group/yzz2/default/legacy/group/projects/vision/rna/')

tss_eRP = c()
dist_eRP = c()

tss_mat = c()
dist_mat = c()
tpm_vec = c()

for (i in c(1:4))
{
print(i)
#i=1
### initiate tss & distal matrix
chrom_id = c(1:19, 'X')
#chrom_id = c(1,2)
### extract state counts
for ( j in chrom_id){
print(j)
Rdat_name = paste('vision_rna_tss2k_ccreunit.chr', toString(j), '.', toString(i), '.1.Rdat', sep='')
load(Rdat_name)
### concatnate ideas state 
tss_mat = rbind(tss_mat, rt$x0)
dist_mat = rbind(dist_mat, rt$nx)
tpm_vec = c(tpm_vec, rt$y)
}
}

### lm regression
tss_model = lm(tpm_vec ~ tss_mat[,-1] - 1)
dist_model = lm(tpm_vec ~ dist_mat[,-1] - 1)
###
tss_eRP = cbind(tss_eRP, tss_model$coefficients)
dist_eRP = cbind(dist_eRP, dist_model$coefficients)
###
tss_eRP = rbind(0, tss_eRP)
rownames(tss_eRP) = c(0:26)
dist_eRP = rbind(0, dist_eRP)
rownames(dist_eRP) = c(0:26)


colnames(tss_eRP) = c('gene_tss0')
colnames(dist_eRP) = c('gene_dist0')

write.table(tss_eRP, 'tss_eRP_2kb.nosplit4G.txt', sep='\t', quote=F, row.names=T, col.names=T)
write.table(dist_eRP, 'dist_eRP_2kb.nosplit4G.txt', sep='\t', quote=F, row.names=T, col.names=T)


### lm regression
tss_model = cor(cbind(tpm_vec, tss_mat))
dist_model = cor(cbind(tpm_vec, dist_mat))
###
tss_eRP = cbind(tss_model[-1,1])
dist_eRP = cbind(dist_model[-1,1])
###
rownames(tss_eRP) = c(0:26)
rownames(dist_eRP) = c(0:26)


colnames(tss_eRP) = c('gene_tss0')
colnames(dist_eRP) = c('gene_dist0')

write.table(tss_eRP, 'tss_eRP_2kb.nosplit4G.cor.txt', sep='\t', quote=F, row.names=T, col.names=T)
write.table(dist_eRP, 'dist_eRP_2kb.nosplit4G.cor.txt', sep='\t', quote=F, row.names=T, col.names=T)


