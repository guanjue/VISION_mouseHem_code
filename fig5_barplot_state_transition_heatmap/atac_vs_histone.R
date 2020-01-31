### get parameters
args = commandArgs(trailingOnly=TRUE)
index_matrix_atac_inputfile = args[1]
index_matrix_histone_inputfile = args[2]
signal_input_list = args[3]
ideas_state_color = args[4]
cREs_IDEASpro_outfile = args[5]


#################################################### 
############ read input files
####################################################
### read signal matrix file
index_matrix_atac_inputfile = 'atac_pk_no0.atac.signal.txt'
index_matrix_histone_inputfile = 'atac_pk_no0.bed.tab.H3K4me3.sort.txt'
signal_input_list = 'signal_list_atac.txt'
signal_input_list_hist = 'signal_list_h3k4me3.txt'
boxplot_name = 'atacpk_histpkporp.h3k4me3.box.pdf'
siglim = 3

print('read signal matrix file')
index_matrix_od = as.matrix(read.table(index_matrix_atac_inputfile, header=FALSE))
index_matrix = index_matrix_od[ , c(5:dim(index_matrix_od)[2]) ]
class(index_matrix) = 'numeric'
pdf('signal_atac.pdf')
hist(index_matrix, breaks=100)
dev.off()

signal_matrix = index_matrix
index_matrix[signal_matrix<=siglim] = 0
index_matrix[signal_matrix>siglim] = 1

histone_matrix_od = as.matrix(read.table(index_matrix_histone_inputfile, header=FALSE))
histone_matrix = histone_matrix_od[ , c(5:dim(histone_matrix_od)[2]) ]
class(histone_matrix) = 'numeric'
pdf('signal_hist.pdf')
hist(histone_matrix, breaks=100)
dev.off()

### sample rows
set.seed(2018)
all_ccRE_num = dim(signal_matrix)[1]
sample_num = 100000
sample_id = sample(all_ccRE_num,sample_num)
signal_matrix_s = signal_matrix[sample_id,]
index_matrix_s = index_matrix[sample_id,]
histone_matrix_s = histone_matrix[sample_id,]
###### read colnames file
print('signal list')
colname_file = read.table(signal_input_list, header=F)
print(colname_file)
colname = colname_file[,2]
colnames(signal_matrix) = colname
colname = colnames(signal_matrix)

colname_file_hist = read.table(signal_input_list_hist, header=F)
print(colname_file_hist)
colname_hist = colname_file_hist[,2]
colnames(histone_matrix) = colname_hist
colname_hist = colnames(histone_matrix)


print(colname)
print(colname_hist)


### cell type index
cpk_hist1_vec = c()
cbg_hist1_vec = c()
t1_pk_hist1_vec = c()
t1_notpk_t2pk_hist1_vec = c()
t1_notpk_hist1_vec = c()

cpk_hist2_vec = c()
cbg_hist2_vec = c()
t2_pk_hist2_vec = c()
t2_notpk_t1pk_hist2_vec = c()
t2_notpk_hist2_vec = c()

ct_pair_vec = c()


k=0

for (i in c("LSK_BM","HPC7","CMP","MEP","CFU_E_ad","CFUMK","MK_imm_ad","GMP")){
	for (j in c("G1E" ,"ER4","ERY_ad","ERY_fl","MONO_BM","NEU","NK_SPL","B_SPL","T_CD4_SPL","T_CD8_SPL")){
		t1col_ct = i
		t2col_ct = j

		if ((sum(colname==t1col_ct)>0) * (sum(colname==t2col_ct)>0) * (sum(colname_hist==t1col_ct)>0) * (sum(colname_hist==t2col_ct)>0) >0){
			k = k+1
			t1_index = index_matrix_s[,colname==t1col_ct]
			t2_index = index_matrix_s[,colname==t2col_ct]

			t1_signal = signal_matrix_s[,colname==t1col_ct]
			t2_signal = signal_matrix_s[,colname==t2col_ct]

			### cell type signal histone
			r1_signal = histone_matrix_s[,colname_hist==t1col_ct]
			r2_signal = histone_matrix_s[,colname_hist==t2col_ct]

			ct_pair = paste(toString(t1col_ct), 'vs', toString(t2col_ct), sep='_')
			print(ct_pair)
			### cpk
			cpk = ((t1_index==1) * (t2_index==1)) == 1
			print(sum(cpk/sample_num))
			cbg = ((t1_index==0) * (t2_index==0)) == 1
			print(sum(cbg/sample_num))
			t1_pk = ((t1_index==1) * (t2_index==0)) == 1
			t1_nopk = (t1_index==0)
			print(sum(t1_pk/sample_num))
			print(sum(t1_nopk/sample_num))
			t2_pk = ((t1_index==0) * (t2_index==1)) == 1
			t2_nopk = (t2_index==0)
			print(sum(t2_pk/sample_num))
			print(sum(t2_nopk/sample_num))

			#cpk_hist1 = (sum((cpk * (r1_signal>siglim))==1)+100) / sum(cpk)+100)
			#cbg_hist1 = (sum((cbg * (r1_signal>siglim))==1)+100) / sum(cbg)+100)
			#t1_pk_hist1 = (sum((t1_pk * (r1_signal>siglim))==1)+100) / sum(t1_pk)+100)
			#t1_notpk_t2pk_hist1 = (sum(((!t1_pk) * t2_pk * (r1_signal>siglim))==1)+100) / sum((!t1_pk) * t2_pk)+100)
			#t1_notpk_hist1 = (sum(((!t1_pk) * (r1_signal>siglim))==1)+100) / sum(!t1_pk)+100)

			#cpk_hist2 = (sum((cpk * (r2_signal>siglim))==1)+100) / sum(cpk)+100)
			#cbg_hist2 = (sum((cbg * (r2_signal>siglim))==1)+100) / sum(cbg)+100)
			#t2_pk_hist2 = (sum((t2_pk * (r2_signal>siglim))==1)+100) / sum(t2_pk)+100)
			#t2_notpk_t1pk_hist2 = (sum(((!t2_pk) * t1_pk * (r2_signal>siglim))==1)+100) / sum((!t2_pk) * t1_pk)+100)
			#t2_notpk_hist2 = (sum(((!t2_pk) * (r2_signal>siglim))==1)+100) / sum(!t2_pk)+100)

			add_num = sample_num*0.01
			cpk_hist1 = (sum((cpk * (r1_signal>siglim))==1)+add_num) / (sum(cpk)*sum((r1_signal>siglim))/sample_num+add_num)
			cbg_hist1 = (sum((cbg * (r1_signal>siglim))==1)+add_num) / (sum(cbg)*sum((r1_signal>siglim))/sample_num+add_num)
			t1_pk_hist1 = (sum((t1_pk * (r1_signal>siglim))==1)+add_num) / (sum(t1_pk)*sum((r1_signal>siglim))/sample_num+add_num)
			t1_notpk_t2pk_hist1 = (sum(((t1_nopk) * t2_pk * (r1_signal>siglim))==1)+add_num) / (sum((t1_nopk)*t2_pk)*sum(r1_signal>siglim)/sample_num+add_num)
			t1_notpk_hist1 = (sum((t1_nopk * (r1_signal>siglim))==1)+add_num) / (sum(t1_nopk)/sample_num*sum((r1_signal>siglim))+add_num)

			cpk_hist2 = (sum((cpk * (r2_signal>siglim))==1)+add_num) / (sum(cpk)*sum((r2_signal>siglim))/sample_num+add_num)
			cbg_hist2 = (sum((cbg * (r2_signal>siglim))==1)+add_num) / (sum(cbg)*sum((r2_signal>siglim))/sample_num+add_num)
			t2_pk_hist2 = (sum((t2_pk * (r2_signal>siglim))==1)+add_num) / (sum(t2_pk)*sum((r2_signal>siglim))/sample_num+add_num)
			t2_notpk_t1pk_hist2 = (sum(((t2_nopk) * t1_pk * (r2_signal>siglim))==1)+add_num) / (sum((t2_nopk)*t1_pk)*sum(r2_signal>siglim)/sample_num+add_num)
			t2_notpk_hist2 = (sum((t2_nopk * (r2_signal>siglim))==1)+add_num) / (sum(t2_nopk)/sample_num*sum((r2_signal>siglim))+add_num)
			
	

			print(cpk_hist1)
			print(cbg_hist1)
			print(paste('t1_pk_hist1:', toString(t1_pk_hist1)))
			print(t1_notpk_t2pk_hist1)
			print(t1_notpk_hist1)

			print(cpk_hist2)
			print(cbg_hist2)
			print(paste('t2_pk_hist2:', toString(t2_pk_hist2)))
			print(t2_notpk_t1pk_hist2)
			print(t2_notpk_hist2)

			ct_pair_vec[k] = ct_pair
			cpk_hist1_vec[k] = cpk_hist1
			cbg_hist1_vec[k] = cbg_hist1
			t1_pk_hist1_vec[k] = t1_pk_hist1
			t1_notpk_t2pk_hist1_vec[k] = t1_notpk_t2pk_hist1
			t1_notpk_hist1_vec[k] = t1_notpk_hist1

			cpk_hist2_vec[k] = cpk_hist2
			cbg_hist2_vec[k] = cbg_hist2
			t2_pk_hist2_vec[k] = t2_pk_hist2
			t2_notpk_t1pk_hist2_vec[k] = t2_notpk_t1pk_hist2
			t2_notpk_hist2_vec[k] = t2_notpk_hist2
		}
	}
}


info_matrix = rbind(cbind(rep('cpk1', length(cpk_hist1_vec)), cpk_hist1_vec),
	cbind(rep('cbg1', length(cbg_hist1_vec)), cbg_hist1_vec),
	cbind(rep('t1pk_t2nopk', length(t1_pk_hist1_vec)), t1_pk_hist1_vec),
	cbind(rep('t1nopk_t2pk', length(t1_notpk_t2pk_hist1_vec)), t1_notpk_t2pk_hist1_vec),
	cbind(rep('cpk2', length(cpk_hist2_vec)), cpk_hist2_vec),
	cbind(rep('cbg2', length(cbg_hist2_vec)), cbg_hist2_vec),
	cbind(rep('t2pk_t1nopk', length(t2_pk_hist2_vec)), t2_pk_hist2_vec),
	cbind(rep('t2nopk_t1pk', length(t2_notpk_t1pk_hist2_vec)), t2_notpk_t1pk_hist2_vec)
)


#png('test.txt')
#boxplot(t1_notpk_hist1_vec)
#dev.off()
#print(table(t1_notpk_hist1_vec))

info_matrix = as.data.frame(info_matrix)
colnames(info_matrix) = c('peaktype', 'sig')

info_matrix$peaktype = factor(info_matrix$peaktype, levels = c('cpk1','cpk2','cbg1','cbg2','t1pk_t2nopk','t2pk_t1nopk','t1nopk_t2pk','t2nopk_t1pk'),ordered = TRUE)
info_matrix[,2] = apply(info_matrix, 1, function(x) as.numeric(x[2]))


#print(info_matrix)


library(ggplot2)

pdf(boxplot_name, width=5, height=2.5)
p = ggplot(data = info_matrix, aes(x=peaktype, y=sig)) 
p = p + geom_boxplot(aes(fill = peaktype))
p = p + geom_point(aes(y=sig, group=peaktype), position = position_dodge(width=0.75))
p = p + scale_fill_manual(values=c("deepskyblue","brown1","deepskyblue","brown1","deepskyblue","brown1",'deepskyblue',"brown1")) #+ theme(panel.background = element_blank(), panel.border = element_rect(colour = "black", fill=NA, size=1.5))
#p = p + stat_compare_means(aes(group = method), label = "p.format", paired = TRUE, method = "t.test")
#p + theme(axis.text.x = element_text(angle = 90, hjust = 1))
p + theme(axis.text.x = element_text(angle = 90, hjust = 1), panel.background = element_blank(), panel.border = element_rect(colour = "black", fill=NA, size=1.5))
plot(p)
dev.off()



siglim_mean = 1

print(mean(r1_signal[(cpk * (r1_signal>siglim_mean))==1]))
print(mean(r1_signal[(cbg * (r1_signal>siglim_mean))==1]))
print(mean(r1_signal[(t1_pk * (r1_signal>siglim_mean))==1]))
print(mean(r2_signal[(t1_pk * (r1_signal>siglim_mean))==1]))

print(mean(r2_signal[(cpk * (r2_signal>siglim_mean))==1]))
print(mean(r2_signal[(cbg * (r2_signal>siglim_mean))==1]))
print(mean(r2_signal[(t2_pk * (r2_signal>siglim_mean))==1]))
print(mean(r1_signal[(t2_pk * (r2_signal>siglim_mean))==1]))

