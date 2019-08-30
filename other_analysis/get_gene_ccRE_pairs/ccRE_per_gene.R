pair_mat = c()

for (i in c(1:4))
{
print(i)
chrom_id = c(1:19, 'X')
#chrom_id = c(1,2)
### extract state counts
for ( j in chrom_id){
print(j)
Rdat_name = paste('vision_rna_tss2k_ccreunit.chr', toString(j), '.', toString(i), '.1.Rdat', sep='')
load(Rdat_name)
### concatnate ideas state 
pair_mat = rbind(pair_mat, rt$pair)
}
}


tss_uniq = unique(pair_mat[,1])

per_gene = dim(pair_mat)[1] / length(tss_uniq)


tss_uniq_02 = unique(pair_mat[pair_mat[,4]>0.2,1])

per_gene_02 = dim(pair_mat[pair_mat[,4]>0.2,])[1] / length(tss_uniq_02)

