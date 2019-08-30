cd /storage/home/gzx103/group/projects/vision/rna


sort -k1,1 -k2,2n rnaTPM.2KB.txt > rnaTPM.2KB.sort.txt
bedtools intersect -a rnaTPM.2KB.sort.txt -b vision_cres.bed -wao > rnaTPM.2KB.ccRE.txt


sort -k1,1 -k2,2n rnaTPM.2MB.txt > rnaTPM.2MB.sort.txt
bedtools intersect -a rnaTPM.2MB.sort.txt -b vision_cres.bed -wao > rnaTPM.2MB.ccRE.txt


R

d=read.table('vision_rna.gene_ccRE.selected.all.reprod_count.WithName.namesort.txt', header=F)

d_02 = d[d[,8]>0.2,]

tss_num = length(table(d_02[,3]))

summary(d)

dim(d_02)[1] / tss_num


