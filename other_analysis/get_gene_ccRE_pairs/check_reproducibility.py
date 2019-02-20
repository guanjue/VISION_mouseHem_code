import numpy as np

################################################################################################
### read 2d array
def read2d_array(filename,dtype_used):
	import numpy as np
	data=open(filename,'r')
	data0=[]
	for records in data:
		tmp = [x.strip() for x in records.split('\t')]
		data0.append(tmp)
	data0 = np.array(data0,dtype=dtype_used)
	data.close()
	return data0

################################################################################################
### write 2d matrix
def write2d_array(array,output):
	r1=open(output,'w')
	for records in array:
		for i in range(0,len(records)-1):
			r1.write(str(records[i])+'\t')
		r1.write(str(records[len(records)-1])+'\n')
	r1.close()


repro_dict = {}
allinfo_dict = {}
uniq_vec = []
for i in range(1,13):
	print(i)
	filename = 'vision_rna.gene_ccRE.c'+str(i)+'.selected.protein_coding.txt'
	d = read2d_array(filename, str)
	for info in d:
		uniq = info[0]+'_'+info[1]+'_'+info[4]
		if not(uniq in repro_dict):
			repro_dict[uniq] = 1
			allinfo_dict[uniq] = info
			uniq_vec.append(uniq)
		else:
			repro_dict[uniq] = repro_dict[uniq]+1


allinfo_mat = np.empty((0, 9))
for uniq in uniq_vec:
	if allinfo_mat.shape[0]%10000 ==0:
		print(allinfo_mat.shape)
	repro_tmp = [repro_dict[uniq]]
	allinfo_tmp = allinfo_dict[uniq]
	allinfo_repro_tmp = np.concatenate((allinfo_tmp,repro_tmp))
	allinfo_mat = np.concatenate((allinfo_mat,allinfo_repro_tmp.reshape(1,9)), axis=0)

result = open('vision_rna.gene_ccRE.selected.protein_coding.reprod_count.txt', 'w')
for uniq in uniq_vec:
	if allinfo_mat.shape[0]%10000 ==0:
		print(allinfo_mat.shape)
	repro_tmp = repro_dict[uniq]
	allinfo_tmp = allinfo_dict[uniq]
	result.write((allinfo_tmp[0])+'\t'+(allinfo_tmp[1])+'\t'+(allinfo_tmp[2])+'\t'+(allinfo_tmp[3])+'\t'+(allinfo_tmp[4])+'\t'+(allinfo_tmp[5])+'\t'+(allinfo_tmp[6])+'\t'+(allinfo_tmp[7])+'\t'+str(repro_tmp)+'\n')

result.close()




allinfo_mat = np.empty((0, 9))
for uniq in uniq_vec:
	if allinfo_mat.shape[0]%10000 ==0:
		print(allinfo_mat.shape)
	repro_tmp = [repro_dict[uniq]]
	allinfo_tmp = allinfo_dict[uniq]
	allinfo_repro_tmp = np.concatenate((allinfo_tmp,repro_tmp))
	allinfo_mat = np.concatenate((allinfo_mat,allinfo_repro_tmp.reshape(1,9)), axis=0)

write2d_array(allinfo_mat, 'vision_rna.gene_ccRE.selected.protein_coding.reprod_count.txt')









