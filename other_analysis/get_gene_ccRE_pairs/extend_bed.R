args = commandArgs(trailingOnly=TRUE)

input_bed = args[1]
expand_win = as.numeric(args[2])
output_file = args[3]

d = read.table(input_bed, header=T, sep=' ')

d_pos = d[,2]
d_pos_exp = cbind(d_pos-expand_win, d_pos+expand_win)
d_pos_exp[d_pos_exp<0]=0

d_sig = d[,5:16]
d_sig_paste = apply(d_sig, 1, function(x) paste(x, collapse='_'))

d_out = cbind(d[,1], d_pos_exp, d[,3:4], d_sig_paste)

d_out = d_out[rowSums(d_pos_exp)>0,]

write.table(d_out, output_file, quote=F, col.names=F, row.names=F, sep='\t')
