for (i in 1:4){
dist1 = read.table(paste('dist.',1,'.sort.bedgraph', sep=''), header=F)
pdf(paste('dist',i,'.hist.pdf', sep=''))
hist(dist1[,4], breaks=50)
dev.off()
}

for (i in 1:4){
dist1 = read.table(paste('tss.',1,'.sort.bedgraph', sep=''), header=F)
pdf(paste('tss',i,'.hist.pdf', sep=''))
hist(dist1[,4], breaks=50)
dev.off()
}


for (i in 1:4){
dist1 = read.table(paste('dist.',1,'.sort.bedgraph', sep=''), header=F)
pdf(paste('dist.no0',i,'.hist.pdf', sep=''))
hist(dist1[,4][dist1[,4]!=0], breaks=100)
dev.off()
}

for (i in 1:4){
dist1 = read.table(paste('tss.',1,'.sort.bedgraph', sep=''), header=F)
pdf(paste('tss.no0',i,'.hist.pdf', sep=''))
hist(dist1[,4][dist1[,4]!=0], breaks=100)
dev.off()
}

