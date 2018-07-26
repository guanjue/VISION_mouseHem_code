### get parameters
args = commandArgs(trailingOnly=TRUE)

intersect_pk_bed_file = args[1]
all_ccRE_pk_bed_file = args[2]
all_TF_pk_bed_file = args[3]
whole_genome_size = as.numeric(args[4])
output_name = args[5]
#2725521370
#2730871774
all_TF_pk_bed = read.table(all_TF_pk_bed_file, header=F, sep='\t')
all_TF_pk_region = sum(all_TF_pk_bed[,3]-all_TF_pk_bed[,2])
#print(all_TF_pk_region)
all_ccRE_pk_bed = read.table(all_ccRE_pk_bed_file, header=F, sep='\t')
all_ccRE_pk_region = sum(all_ccRE_pk_bed[,3]-all_ccRE_pk_bed[,2])
print(all_ccRE_pk_region)
intersect_pk_bed = read.table(intersect_pk_bed_file, header=F, sep='\t')
intersect_pk_region = sum(intersect_pk_bed[,3]-intersect_pk_bed[,2])
#print(intersect_pk_region)
#enrichment = (intersect_pk_region+100.0) / ((all_TF_pk_region+0.0) * (all_ccRE_pk_region+0.0) / whole_genome_size + 100.0)
enrichment = (dim(intersect_pk_bed)[1]+10.0) / ((dim(all_TF_pk_bed)[1]+0.0) * (all_ccRE_pk_region+0.0) / whole_genome_size + 10.0)
print('result')
print(dim(intersect_pk_bed)[1])
print((dim(all_TF_pk_bed)[1]+0.0) * (all_ccRE_pk_region+0.0) / whole_genome_size)
print(enrichment)
write.table(c(all_ccRE_pk_bed_file, enrichment), paste(output_name, '.enrichment.txt', sep=''), quote=F, sep='\t', col.names=F, row.names=F)

