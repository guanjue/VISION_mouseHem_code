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
d0 = open('all_converted.sort.txt', 'r')
d1 = []

### read all mat
for records in d0:
	tmp = [x.strip() for x in records.split('\t')]
	d1.append(tmp)
d1 = np.array(d1,dtype=str)
d0.close()

### read gene group
def get_gene_group_dict(input_name):
	### read input
	g0 = open(input_name, 'r')
	g0a = []
	for records in g0:
		tmp = [x.strip() for x in records.split('\t')]
		g0a.append(tmp)
	g0a = np.array(g0a,dtype=str)
	g0.close()
	### get dict
	g0_dict = {}
	for records in g0a:
		if not (records[0] in g0_dict):
			g0_dict[records[0]] = '1'
	### output
	return g0_dict


g1_dict = get_gene_group_dict('gene.group.1.txt')
g2_dict = get_gene_group_dict('gene.group.2.txt')
g3_dict = get_gene_group_dict('gene.group.3.txt')
g4_dict = get_gene_group_dict('gene.group.4.txt')


d1_g = []

for records in d1:
	### get tss
	tss = records[0]+'_'+records[1]
	if tss in g1_dict:
		g = '1'
		d1_g.append([records[0], records[1], records[2], records[3], records[4], records[5], records[6], g])
	elif tss in g2_dict:
		g = '2'
		d1_g.append([records[0], records[1], records[2], records[3], records[4], records[5], records[6], g])
	elif tss in g3_dict:
		g = '3'
		d1_g.append([records[0], records[1], records[2], records[3], records[4], records[5], records[6], g])
	elif tss in g4_dict:
		g = '4'
		d1_g.append([records[0], records[1], records[2], records[3], records[4], records[5], records[6], g])

d1_g = np.array(d1_g)
write2d_array(d1_g, 'all_converted.sort.g.txt')

