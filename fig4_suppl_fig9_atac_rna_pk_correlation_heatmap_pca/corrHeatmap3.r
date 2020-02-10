args <- commandArgs(trailingOnly = TRUE)
intab <- args[1]

#normalize data first => filter all zero, add min, log2, quantile?
data=read.table(intab, sep="\t", header=TRUE, row.names=1)

c=cor(data, method="spearman") #pearson is default, use spearman to match ENCODE
data=as.matrix(data)
data=as.numeric(data)

cr=round(c, digits=2)
library("gplots")
library("RColorBrewer")
my_palette <- c("#1000FFFF","#0047FFFF","#00C8FFFF","#00FFE0FF","#00FF89FF","#00FF33FF","#99FF00FF","#AFFF00FF","#FFD800FF","#FFB800FF","#FFAD00FF","#FF8100FF","#FF5600FF","#FF2B00FF","#FF0000FF")

#breaks for RNA
breaks = c(seq(.7,1,.02))
dist.pear <- function(x) as.dist(1-x)
heatmap.2(c, dendrogram="row", col=my_palette, breaks=breaks, trace="none", revC=TRUE, symm=TRUE, distfun=dist.pear, symkey=FALSE, keysize=1.3, key.title=NA, key.par=list(cex=1)) 

q()
