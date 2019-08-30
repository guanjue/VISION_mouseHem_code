d = read.table('tad.net.txt', header = F)
c=d[,2]
pdf('edge_count_hist.pdf')
hist(d[,2], breaks=50)
dev.off()
ct = quantile(c, 0.95)
dt = d[c>ct,]
write.table(dt, 'tad.net.sub.txt', sep='\t', quote=F, col.names=F, row.names=F)


d = read.table('tad.netenrich.txt', header = F)
c=d[,2]
pdf('edge_enrich_count_hist.pdf')
hist(d[,2], breaks=50)
dev.off()
ct = quantile(c, 0.95)
dt = d[c>ct,]
write.table(dt, 'tad.netenrich.sub.txt', sep='\t', quote=F, col.names=F, row.names=F)


