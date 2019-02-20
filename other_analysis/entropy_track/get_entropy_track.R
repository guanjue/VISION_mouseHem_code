
d = read.table('pknorm_2_16lim_ref1mo_0424_lesshet.state', header=F, sep=' ')

count_to_entropy = function(x){
	c_tmp = table(x)
	return(-sum(c_tmp * log(c_tmp)))
}

d_bed = d[,2:4]

d_label = d[,-c(1:4, dim(d)[2])]

d_label_entropy = apply(d_label, 1, count_to_entropy)


d_label_entropy_bedgraph = cbind(d_bed, d_label_entropy)

write.table(d_label_entropy_bedgraph, 'pknorm_2_16lim_ref1mo_0424_lesshet.entropy.bedgraph', quote=F, sep='\t', col.names=F, row.names=F)


d_label_entropy_bedgraph = cbind(d_bed, -d_label_entropy)

write.table(d_label_entropy_bedgraph, 'pknorm_2_16lim_ref1mo_0424_lesshet.neg_entropy.bedgraph', quote=F, sep='\t', col.names=F, row.names=F)




