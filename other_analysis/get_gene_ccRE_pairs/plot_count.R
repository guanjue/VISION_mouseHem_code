### read count table
d1 = read.table('all_converted.sort.g1.uniq.count.txt', header=F)
d2 = read.table('all_converted.sort.g2.uniq.count.txt', header=F)
d3 = read.table('all_converted.sort.g3.uniq.count.txt', header=F)
d4 = read.table('all_converted.sort.g4.uniq.count.txt', header=F)

### get counts
d1_allccRE = d1[,5]
d1_cor02ccRE = d1[,6]
d1_selectccRE = d1[,7]

d2_allccRE = d2[,5]
d2_cor02ccRE = d2[,6]
d2_selectccRE = d2[,7]

d3_allccRE = d3[,5]
d3_cor02ccRE = d3[,6]
d3_selectccRE = d3[,7]

d4_allccRE = d4[,5]
d4_cor02ccRE = d4[,6]
d4_selectccRE = d4[,7]

### plot boxplot

pdf('ccREcount.boxplot.pdf', width=28)
par(mfrow=c(1,4))
boxplot(cbind(d1_allccRE, d1_cor02ccRE, d1_selectccRE), ylim = c(0, 500), col=(c('white','blue','red')), main = paste(toString(round(mean(d1_allccRE), 1)), toString(round(mean(d1_cor02ccRE), 1)), toString(round(mean(d1_selectccRE), 1)), sep='---' ) )
boxplot(cbind(d2_allccRE, d2_cor02ccRE, d2_selectccRE), ylim = c(0, 500), col=(c('white','blue','red')), main = paste(toString(round(mean(d1_allccRE), 1)), toString(round(mean(d2_cor02ccRE), 1)), toString(round(mean(d2_selectccRE), 1)), sep='---' ) )
boxplot(cbind(d3_allccRE, d3_cor02ccRE, d3_selectccRE), ylim = c(0, 500), col=(c('white','blue','red')), main = paste(toString(round(mean(d1_allccRE), 1)), toString(round(mean(d3_cor02ccRE), 1)), toString(round(mean(d3_selectccRE), 1)), sep='---' ) )
boxplot(cbind(d4_allccRE, d4_cor02ccRE, d4_selectccRE), ylim = c(0, 500), col=(c('white','blue','red')), main = paste(toString(round(mean(d1_allccRE), 1)), toString(round(mean(d4_cor02ccRE), 1)), toString(round(mean(d4_selectccRE), 1)), sep='---' ) )
dev.off()



