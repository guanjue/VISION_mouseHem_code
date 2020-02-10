#setwd('/gpfs/group/yzz2/default/legacy/group/projects/vision/rna/')
### get parameters
args = commandArgs(trailingOnly=TRUE)

wddir = args[1]
setwd(wddir)

tss_coefficient = c()
dist_coefficient = c()

for (i in c(1:4))
{
print(i)
#i=1
### initiate tss & distal matrix
tss_mat = c()
dist_mat = c()
tpm_vec = c()
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
### lm regression
tss_model = lm(tpm_vec ~ tss_mat - 1)
dist_model = lm(tpm_vec ~ dist_mat - 1)
###
tss_coefficient = cbind(tss_coefficient, tss_model$coefficients)
dist_coefficient = cbind(dist_coefficient, dist_model$coefficients)
}

#tss_coefficient = rbind(rep(0.0,4), tss_coefficient)
#dist_coefficient = rbind(rep(0.0,4), dist_coefficient)

colnames(tss_coefficient) = c('gene_c1', 'gene_c2', 'gene_c3', 'gene_c4')
colnames(dist_coefficient) = c('gene_c1', 'gene_c2', 'gene_c3', 'gene_c4')

rownames(tss_coefficient) = c(0:26)
rownames(dist_coefficient) = c(0:26)

write.table(tss_coefficient, 'tss_eRP_2kb.with0.txt', sep='\t', quote=F, row.names=T, col.names=T)
write.table(dist_coefficient, 'dist_eRP_2kb.with0.txt', sep='\t', quote=F, row.names=T, col.names=T)



