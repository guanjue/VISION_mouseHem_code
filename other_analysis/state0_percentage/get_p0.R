setwd('~/group/projects/vision/pknorm_2_16lim_ref1mo_sample1mo_0424')

d=read.table('pknorm_2_16lim_ref1mo_0424_lesshet.state')


coln = c('ID', 'CHR', 'POSst', 'POSed', 'B_SPL', 'CFU_E_ad', 'CFUMK', 'CLP', 'CMP', 'ER4', 'ERY_ad', 'ERY_fl', 'G1E', 'GMP', 'HPC7', 'LSK_BM', 'MEP', 'MK_imm_ad', 'MK_mat_fl', 'MONO_BM', 'NEU', 'NK_SPL', 'T_CD4_SPL', 'T_CD8_SPL', 'PosClass')

bin_num = dim(d)[1]

p0 = c()
for (i in c(1:20)){
	print(coln[i+4])
	p0_tmp = sum(d[,i+4]==0)/bin_num
	print(p0_tmp)
	p0 = c(p0, p0_tmp)
}

p0_all = cbind(coln[-c(1:4,25)], p0)
print(p0_all[order(p0),])
write.table(p0_all, 'celltype_s0_p.txt', sep=' ', quote=F)


all_0 = apply(d,1,function(x) sum(as.numeric(x[5:24])) )
all_0_p = sum(all_0==0)/bin_num
print('all regions in state 0: ')
print(all_0_p)

write.table(all_0_p, 'allcelltype_s0_p.txt', sep=' ', quote=F)
