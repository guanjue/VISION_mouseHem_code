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
def get_uniq(input_name):
	d0 = open(input_name+'.txt', 'r')
	d1 = []
	print(input_name)
	### read all mat
	for records in d0:
		tmp = [x.strip() for x in records.split('\t')]
		d1.append(tmp)
	d1 = np.array(d1,dtype=str)
	d0.close()

	### 
	d1_uniq = []
	d1_uniq_dict = {}

	tss_ccRE_cor = {}
	tss_ccRE_repr = {}

	for records in d1:
		### get tss
		tss_ccRE = records[0]+'_'+records[1]+'_'+records[2]+'_'+records[3]+'_'+records[4]+'_'+records[7]
		if not (tss_ccRE in d1_uniq_dict):
			d1_uniq_dict[tss_ccRE] = ''
			d1_uniq.append(tss_ccRE)
			### add first info
			tss_ccRE_cor[tss_ccRE] = records[5]
			tss_ccRE_repr[tss_ccRE] = records[6]
		elif (tss_ccRE_cor[tss_ccRE] == '.') & (records[5] == '.'):
			print('skip')
		elif (tss_ccRE_cor[tss_ccRE] != '.') & (records[5] == '.'):
			print('skip')
		elif (tss_ccRE_cor[tss_ccRE] == '.') & (records[5] != '.'):
			tss_ccRE_cor[tss_ccRE] = records[5]
			tss_ccRE_repr[tss_ccRE] = records[6]
		elif (tss_ccRE_cor[tss_ccRE] != '.') & (records[5] != '.'):
			### check repr
			if (tss_ccRE_repr[tss_ccRE]=='0') & (records[6]=='1'):
				tss_ccRE_cor[tss_ccRE] = records[5]
				tss_ccRE_repr[tss_ccRE] = records[6]			
			elif (tss_ccRE_repr[tss_ccRE]=='1') & (records[6]=='0'):
				print('skip')
			else:
				### check cor
				if float(tss_ccRE_cor[tss_ccRE]) < float(records[5]):
					tss_ccRE_cor[tss_ccRE] = records[5]
					tss_ccRE_repr[tss_ccRE] = records[6]


	d_new = []		
	for tss_ccRE in d1_uniq:
		tss_ccRE_vec = tss_ccRE.split('_')
		tmp_cor = tss_ccRE_cor[tss_ccRE]
		tmp_repr = tss_ccRE_repr[tss_ccRE]
		d_new.append([tss_ccRE_vec[0], tss_ccRE_vec[1], tss_ccRE_vec[2], tss_ccRE_vec[3], tss_ccRE_vec[4], tss_ccRE_vec[5], tmp_cor, tmp_repr])

	d_new = np.array(d_new)
	write2d_array(d_new, input_name+'.uniq.txt')

get_uniq('all_converted.sort.g1')
get_uniq('all_converted.sort.g2')
get_uniq('all_converted.sort.g3')
get_uniq('all_converted.sort.g4')




