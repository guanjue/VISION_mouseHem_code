args = commandArgs(trailingOnly=TRUE)
ct = args[1]
tss_cof = args[2]
dist_cof = args[3]

d = read.table(paste('ccRE.',ct,'.state.2.all.statepercent.txt', sep=''), header=F)
dpk = read.table(paste('ccRE.',ct,'.state.2.all.bed', sep=''), header=F)

tss_eRPs = read.table(tss_cof, header=T)
tss_eRPs = t(apply(tss_eRPs,1,function(x) x-as.matrix(tss_eRPs[1,])))

dist_eRPs = read.table(dist_cof, header=T)
dist_eRPs = t(apply(dist_eRPs,1,function(x) x-as.matrix(dist_eRPs[1,])))


tss_eRP_mat = c()
for (i in 1:4){
	print(i)
tss_eRP_1 = (apply(d, 1, function(x) sum(x/sum(x)*tss_eRPs[,i])))
tss_eRP_mat = cbind(tss_eRP_mat, tss_eRP_1)
}

dist_eRP_mat = c()
for (i in 1:4){
	print(i)
dist_eRP_1 = (apply(d, 1, function(x) sum(x/sum(x)*dist_eRPs[,i])))
dist_eRP_mat = cbind(dist_eRP_mat, dist_eRP_1)
}

all_eRPs = cbind(tss_eRP_mat, dist_eRP_mat)
output_file = c('tss.1.bedgraph','tss.2.bedgraph','tss.3.bedgraph','tss.4.bedgraph',
	'dist.1.bedgraph','dist.2.bedgraph','dist.3.bedgraph','dist.4.bedgraph')

output_file = cbind(rep(ct,8), c('tss.1.bedgraph','tss.2.bedgraph','tss.3.bedgraph','tss.4.bedgraph',
	'dist.1.bedgraph','dist.2.bedgraph','dist.3.bedgraph','dist.4.bedgraph'))

output_file = apply(output_file, 1, function(x) paste(x[1], x[2], sep='.'))


for (i in 1:dim(all_eRPs)[2]){
	bedgraph = cbind(dpk[,1:3], all_eRPs[,i])
	write.table(bedgraph, output_file[i], quote=F, sep='\t', col.names=F, row.names=F)
}


