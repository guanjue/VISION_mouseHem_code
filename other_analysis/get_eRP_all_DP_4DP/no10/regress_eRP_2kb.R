setwd('/gpfs/group/yzz2/default/legacy/group/projects/vision/rna/')

ccRE_eRP = c()
ccRE_mat = c()
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
Rdat_name = paste('vision_rna_tss2k_ccreunit_folder_used0509/vision_rna_tss2k_ccreunit.chr', toString(j), '.', toString(i), '.1.Rdat', sep='')
load(Rdat_name)
### concatnate ideas state 
tss_mat = rt$x0 * 2000
dist_mat = rt$nx * (2000000-2000)
dist_mat_all = (tss_mat + dist_mat) / 2000000
ccRE_mat = rbind(ccRE_mat, dist_mat_all)
tpm_vec = c(tpm_vec, rt$y)
}
}

### lm regression
ccRE_model = lm(tpm_vec ~ ccRE_mat[,-1] - 1)
###
ccRE_eRP = cbind(c(0.0, ccRE_model$coefficients))
###
#ccRE_eRP = rbind(0, ccRE_eRP)
rownames(ccRE_eRP) = c(0:26)

colnames(ccRE_eRP) = c('gene_ccRE0')

write.table(ccRE_eRP, 'eRP_all.txt', sep='\t', quote=F, row.names=T, col.names=T)

