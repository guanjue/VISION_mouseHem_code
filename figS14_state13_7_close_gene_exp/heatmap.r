args <- commandArgs(trailingOnly = TRUE)
intab <- args[1]
intab2 <- args[2]

data=read.table(pipe(paste("cut -f 4- ", intab, sep="")), sep="\t", header=TRUE, row.names=1)
data2=read.table(pipe(paste("cut -f 4- ", intab2, sep="")), sep="\t", header=TRUE, row.names=1)

#x=cells y=genes, group genes by class based on state
data=data.matrix(data)
data2=data.matrix(data2)

library("gplots")
library("RColorBrewer")

pdf("heatmap13_13.pdf")
breaks = c(seq(-10,10,.2))
heatmap.2(data, hclustfun = function(x) hclust(x, method="ward.D"), dendrogram="none", trace="none", col=colorRampPalette(c("white","white","red")), Colv=FALSE, symkey=FALSE, keysize=1.3, key.title=NA, key.par=list(cex=1), labRow=FALSE, breaks=breaks)
dev.off()
pdf("heatmap13_7.pdf")
heatmap.2(data2, hclustfun = function(x) hclust(x, method="ward.D"), dendrogram="none", trace="none", col=colorRampPalette(c("white","white","red")), Colv=FALSE, symkey=FALSE, keysize=1.3, key.title=NA, key.par=list(cex=1), labRow=FALSE, breaks=breaks)
dev.off()

q()

