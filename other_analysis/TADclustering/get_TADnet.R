
net_all = c()
nodes_all = c()
pk_all = c()

for (i in c(1:856)){
print(i)
filename = paste('module_pk_list/TAD_and_bd.', i, '.state.bed', sep='')
d = read.table(filename, header=F, sep='\t')
for (j in c(1:(dim(d)[1]-1))){
edge = paste(d[j,4], ',', d[j+1,4], sep='')
net_all = rbind(net_all, edge)
nodes_all = rbind(nodes_all, as.character(d[j,4]))
pk_all = rbind(pk_all, d[j,])
}
}


bd3_tad6_pk = pk_all[net_all=="bd_3,tad_6",]


net_all_counts = table(net_all)
write.table(cbind(net_all_counts), 'tad.net.txt', quote=F, sep='\t', col.names=F)

node_all_counts = table(nodes_all)
write.table(cbind(node_all_counts), 'tad.node.txt', quote=F, sep='\t', col.names=F)


edge_name = rownames(cbind(net_all_counts))
edge_enrich = c()
for (i in 1:length(net_all_counts)){
edge1 = unlist(strsplit(edge_name[i], ','))
edge_c = net_all_counts[i]
exp_c = node_all_counts[rownames(node_all_counts)==edge1[1]] * node_all_counts[rownames(node_all_counts)==edge1[2]] / sum(node_all_counts)
edge_enrich_tmp = (edge_c+10) / (exp_c+10)
edge_enrich = rbind(edge_enrich, c(edge_name[i], edge_enrich_tmp))
}

write.table(edge_enrich, 'tad.netenrich.txt', quote=F, sep='\t', col.names=F, row.names=F)



pdf('edge_count.pdf', width=14)
par(mfrow=c(1,2))
hist(net_all_counts, breaks=50)
hist(log2(net_all_counts), breaks=50)
dev.off()


