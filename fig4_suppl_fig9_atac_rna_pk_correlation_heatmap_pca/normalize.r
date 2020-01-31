args <- commandArgs(trailingOnly = TRUE)
intab <- args[1]
outtab <- args[2]
data <- read.table(intab, sep="\t", header=TRUE, row.names=1)
#remove transcript IDs if present
if(any(grepl("transcript", colnames(data)))) {
   data=data[, -grep("transcript", colnames(data))]
}
a=min(data[data!=0])
dim(data)
#filter genes with all zero
data=data[rowSums(data[,-1])>0,]
dim(data)
#take log2 of value+min 
ldata=log2(data + a)
#quantile
library(limma)
limmaRes = normalizeQuantiles(ldata)
limmaRes = round(limmaRes, digits=3)
write.table(limmaRes, file=outtab, row.names = TRUE ,col.names = TRUE,quote = FALSE, sep="\t")

