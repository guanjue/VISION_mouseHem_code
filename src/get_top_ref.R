library(metap)

### get parameters
args = commandArgs(trailingOnly=TRUE)

tail = args[1]
input_folder = args[2]
output_name = args[3]

### extract filenames of the cell marker
ref_file_list = list.files(input_folder, pattern=paste(tail, '$', sep='') )
print(ref_file_list)

file_list = c()
FRiP_list = c()
SNR_list = c()
for (i in c(1:length(ref_file_list))){
	print(i)
	print(ref_file_list[i])
	file=ref_file_list[i]
	parameters = read.table(file, header=F)
	### FRiP score is in the 5th column
	file_list[i] = as.character(parameters[1,])
	FRiP_list[i] = as.numeric(as.character(parameters[2,]))
	SNR_list[i] = as.numeric(as.character(parameters[3,]))
}


print(FRiP_list)
print(which.max(SNR_list))
print(SNR_list)
print(which.max(SNR_list))


ref_file = file_list[which.max(SNR_list)]
frip_ref = FRiP_list[which.max(SNR_list)]
SNR_ref = SNR_list[which.max(SNR_list)]

write.table(c(ref_file, frip_ref, SNR_ref), output_name, sep='\t', quote=F, col.names=F, row.names=F)

pknorm_list = cbind(rep(ref_file, length(file_list)), file_list)
write.table(pknorm_list, paste(output_name, '.info.txt', sep=''), sep='\t', quote=F, col.names=F, row.names=F)

