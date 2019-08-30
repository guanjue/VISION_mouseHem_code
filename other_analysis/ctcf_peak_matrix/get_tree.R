ct = read.table('ct.list.txt')

d=read.table('chr11_loci.idsort.mat.S3norm_NBP.txt', header=F)
dsig = d[,-c(1:4)]
colnames(dsig) = ct[,1]
hc = hclust(dist(t(dsig)))
pdf('chr11_loci.idsort.mat.S3norm_NBP.tree.pdf')
plot(hc)
dev.off()

d=read.table('chr11_loci.idsort.mat.NBP.txt', header=F)
dsig = d[,-c(1:4)]
colnames(dsig) = ct[,1]
hc = hclust(dist(t(dsig)))
pdf('chr11_loci.idsort.mat.NBP.tree.pdf')
plot(hc)
dev.off()


d=read.table('chr11_loci.idsort.mat.S3norm_rc.txt', header=F)
dsig = d[,-c(1:4)]
colnames(dsig) = ct[,1]
hc = hclust(dist(t(dsig)))
pdf('chr11_loci.idsort.mat.S3norm_rc.tree.pdf')
plot(hc)
dev.off()


