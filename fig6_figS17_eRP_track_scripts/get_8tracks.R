args = commandArgs(trailingOnly=TRUE)
input_file = args[1]
output_file = args[2]


d = read.table(input_file, header=F)
pks_list_all = d[,4]
pks_list = unique(pks_list_all)

### get pk_len
pk_len = c()
for (i in 1:length(pks_list)){
	if ((i%%1000)==0){
		print(i)
	}
	pk_info = unlist(strsplit(as.character(pks_list[i]), '_'))
	pk_len[i] = as.numeric(pk_info[3])-as.numeric(pk_info[2])
}


### 
len_state_mat = c()

for (i in 1:length(pks_list)){
if ((i%%10)==0){
	print(i)
}
len_state_vec = rep(0, 27)
all_len_tmp = pk_len[i]
pk_list_tmp = d[pks_list_all==pks_list[i],c(2,3,5)]
#print(pk_list_tmp[,3])
#print(as.character(pk_list_tmp[,3]))
#print(as.numeric(as.character(pk_list_tmp[,3])))
pk_list_tmp_state = as.numeric(as.character(pk_list_tmp[,3]))
pk_list_tmp_dist = pk_list_tmp[,2]-pk_list_tmp[,1]
if (dim(pk_list_tmp)[1]>0){
for (j in 1:dim(pk_list_tmp)[1]){
len_state_vec[pk_list_tmp_state[j]+1] = len_state_vec[pk_list_tmp_state[j]+1] + pk_list_tmp_dist[j]
}
}
len_state_mat = rbind(len_state_mat, len_state_vec)
}


write.table(len_state_mat, output_file, quote=F, sep='\t', col.names=F, row.names=F)



#tss_eRPs = read.table('~/group/projects/vision/rna/tss_eRP_2kb.with0.txt', header=T)
#dist_eRPs = read.table('~/group/projects/vision/rna/dist_eRP_2kb.with0.txt', header=T)





