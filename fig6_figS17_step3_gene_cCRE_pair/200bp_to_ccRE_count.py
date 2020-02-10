import os
import os.path
import numpy as np
from subprocess import call
from collections import Counter
import decimal


################################################################################################
### read 2d array
def read2d_array(filename,dtype_used):
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

################################################################################################
###
def count_ccRE_pre_gene(input_name):
	d0 = read2d_array(input_name+'.uniq.txt', str)
	d0_dict_all_ccRE = {}
	d0_dict_cor02_ccRE = {}
	d0_dict_select_ccRE = {}
	d0_vec = []

	for records in d0:
		gene = records[2]
		if not (gene in d0_dict_all_ccRE):
			d0_vec.append(records)
			d0_dict_all_ccRE[gene] = 1
			### get cor
			cor_tmp = records[6]
			repr_tmp = records[7]
			### check correlation
			if cor_tmp!='.':
				d0_dict_cor02_ccRE[gene] = 1
				### check selection
				if (repr_tmp != '0') and (repr_tmp != '1'):
					d0_dict_select_ccRE[gene] = 1
				else:
					d0_dict_select_ccRE[gene] = 0
			else:
				d0_dict_cor02_ccRE[gene] = 0
				d0_dict_select_ccRE[gene] = 0
		else:
			d0_dict_all_ccRE[gene] = d0_dict_all_ccRE[gene]+1
			### get cor
			cor_tmp = records[6]
			repr_tmp = records[7]
			### check correlation
			if cor_tmp!='.':
				d0_dict_cor02_ccRE[gene] = d0_dict_cor02_ccRE[gene]+1
				### check selection
				if (repr_tmp != '0') and (repr_tmp != '1'):
					d0_dict_select_ccRE[gene] = d0_dict_select_ccRE[gene]+1
				else:
					d0_dict_select_ccRE[gene] = d0_dict_select_ccRE[gene]+0
			else:
				d0_dict_cor02_ccRE[gene] = d0_dict_cor02_ccRE[gene]+0
				d0_dict_select_ccRE[gene] = d0_dict_select_ccRE[gene]+0

	### output
	dnew = []
	for records in d0_vec:
		gene = records[2]
		all_count_tmp = d0_dict_all_ccRE[gene]
		cor02_count_tmp = d0_dict_cor02_ccRE[gene]
		select_count_tmp = d0_dict_select_ccRE[gene]
		dnew.append([records[0], records[1], records[2], records[5], all_count_tmp, cor02_count_tmp, select_count_tmp])
	dnew = np.array(dnew)
	write2d_array(dnew, input_name+'.uniq.count.txt')


count_ccRE_pre_gene('all_converted.sort.g1')
count_ccRE_pre_gene('all_converted.sort.g2')
count_ccRE_pre_gene('all_converted.sort.g3')
count_ccRE_pre_gene('all_converted.sort.g4')







