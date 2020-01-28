### get parameters
args = commandArgs(trailingOnly=TRUE)

input_HTcount = args[1]
output_HTcount = args[2]
output_TPM = args[3]

#d=read.table('rnaHtseqCounts_withcoordinates.txt', header=T)
d=read.table(input_HTcount, header=T)

dsig = d[,-c(1:7)]


#write.table(cbind(d[,1:6],dsig),'rnaHtseqCounts_withcoordinates.0.txt', quote=F, col.names=T, row.names=F, sep='\t')
write.table(cbind(d[,1:6],dsig),output_HTcount, quote=F, col.names=T, row.names=F, sep='\t')

d_len = d[,3]-d[,2]
dsig_transcripts = dsig/d_len

dsig_transcripts_colsum = colSums(dsig_transcripts)

dsig_tpm = t(apply(dsig_transcripts, 1, function(x) x/dsig_transcripts_colsum*1000000))

dsig_tpm_log2 = log2(dsig_tpm+1e-3)

#write.table(cbind(d[,1:6],dsig),'rnaTPM_withcoordinates.0.txt', quote=F, col.names=T, row.names=F, sep='\t')
write.table(cbind(d[,1:6],dsig_tpm_log2),output_TPM, quote=F, col.names=T, row.names=F, sep='\t')


