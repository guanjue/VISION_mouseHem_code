import os
import os.path
import numpy as np
from subprocess import call
from collections import Counter
import decimal

### get tss vs gene_nameid matching table
d0 = open('gencode.vM4.tss.tss2gene0.txt', 'r')
d0_dict = {}

for records in d0:
	tmp = [x.strip() for x in records.split('\t')]
	if not (tmp[0] in d0_dict):
		d0_dict[tmp[0]] = tmp[1]
d0.close()

### 
d1 = open('vision_rna.gene_ccRE.selected.all.reprod_count.txt','r')
d1a = []
d1new = []

for records in d1:
	tmp = [x.strip() for x in records.split('\t')]
	#tss = tmp[0]+'_'+str(int(float(tmp[1])))
	tss = tmp[0]+'_'+tmp[1]
	#print(tss)
	gene = d0_dict[tss]
	d1new.append([tmp[0], tmp[1], gene, tmp[2], tmp[3], tmp[4], tmp[5], tmp[6], tmp[7], tmp[8]])

d1.close()
d1new = np.array(d1new)

################################################################################################
### write 2d matrix
def write2d_array(array,output):
	r1=open(output,'w')
	for records in array:
		for i in range(0,len(records)-1):
			r1.write(str(records[i])+'\t')
		r1.write(str(records[len(records)-1])+'\n')
	r1.close()

write2d_array(d1new, 'vision_rna.gene_ccRE.selected.all.reprod_count.WithName1209.txt')


