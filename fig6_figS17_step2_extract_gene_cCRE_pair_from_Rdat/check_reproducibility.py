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
	filename = 'vision_rna.gene_ccRE.c'+str(i)+'.selected.all.txt'
	d = read2d_array(filename, str)
	for info in d:
		uniq = info[0]+'_'+info[1]+'_'+info[4]
		if not(uniq in repro_dict):
			repro_dict[uniq] = ''
			allinfo_dict[uniq] = info
			uniq_vec.append(uniq)

for i in range(1,13):
	print(i)
	filename = 'vision_rna.gene_ccRE.c'+str(i)+'.selected.all.txt'
	d = read2d_array(filename, str)
	repro_dict_tmp = {}
	for info in d:
		uniq = info[0]+'_'+info[1]+'_'+info[4]
		if not(uniq in repro_dict_tmp):
			repro_dict_tmp[uniq] = info[7]
	for uniq0 in uniq_vec:
		if uniq0 in repro_dict_tmp:
			repro_dict[uniq0] = repro_dict[uniq0] + repro_dict_tmp[uniq0]
		else:
			repro_dict[uniq0] = repro_dict[uniq0] + '0'


result = open('vision_rna.gene_ccRE.selected.all.reprod_count.txt', 'w')
for uniq in uniq_vec:
	repro_tmp = repro_dict[uniq]
	allinfo_tmp = allinfo_dict[uniq]
	result.write((allinfo_tmp[0])+'\t'+(allinfo_tmp[1])+'\t'+(allinfo_tmp[2])+'\t'+(allinfo_tmp[3])+'\t'+(allinfo_tmp[4])+'\t'+(allinfo_tmp[5])+'\t'+(allinfo_tmp[6])+'\t'+(allinfo_tmp[7])+'\t'+str(repro_tmp)+'\n')

result.close()


