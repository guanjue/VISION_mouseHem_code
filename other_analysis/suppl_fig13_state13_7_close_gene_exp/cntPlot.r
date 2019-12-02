args <- commandArgs(trailingOnly = TRUE)
intab <- args[1]

data=read.table(intab, sep="\t", header=TRUE)

summary(data)
order=c("X50k_13", "X50k_7", "X100k_13", "X100k_7")
t.test(data[,1], data[,3])
t.test(data[,2], data[,4])

library(ggplot2)
pdf("geneCnts.pdf")
p=ggplot(data=data) + geom_violin(aes(x="X50k_13", y=data[,1])) + geom_violin(aes(x="X100k_13", y=data[,2])) + geom_violin(aes(x="X50k_7", y=data[,3])) + geom_violin(aes(x="X100k_7", y=data[,4])) + scale_x_discrete(limits=order) + labs (y = "Count of genes near cCRE")
p
dev.off()

q()

