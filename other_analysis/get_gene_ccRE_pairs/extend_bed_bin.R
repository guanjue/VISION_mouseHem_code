args = commandArgs(trailingOnly=TRUE)

input_bed = args[1]
expand_win1 = as.numeric(args[2])
expand_win2 = as.numeric(args[3])
output_file = args[4]

d = read.table(input_bed, header=T, sep=' ')

d_pos = d[,2]
d_pos_exp1 = cbind(d_pos+expand_win1, d_pos+expand_win2)
d_pos_exp1[d_pos_exp1<0]=0
d_pos_exp2 = cbind(d_pos-expand_win2, d_pos-expand_win1)
d_pos_exp2[d_pos_exp2<0]=0

d_sig = d[,5:16]
d_sig_paste = apply(d_sig, 1, function(x) paste(x, collapse='_'))

d_out1 = cbind(d[,1], d_pos_exp1, d[,3:4], d_sig_paste)
d_out1 = d_out1[rowSums(d_pos_exp1)>0,]
d_out2 = cbind(d[,1], d_pos_exp2, d[,3:4], d_sig_paste)
d_out2 = d_out2[rowSums(d_pos_exp2)>0,]

d_out = rbind(d_out1, d_out2)

write.table(d_out, output_file, quote=F, col.names=F, row.names=F, sep='\t')
