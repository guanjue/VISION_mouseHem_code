import numpy as np

#ccRE0 = open('mouse_allChr_wTAD_022519.bed', 'r')
ccRE0 = open('mouse_allChr_wTAD_022519.repeat.bed', 'r')
ccRE1 = []

for records in ccRE0:
	tmp = [x.strip() for x in records.split('\t')]
	ccRE1.append(tmp[3])


ccRE_eRP0 = open('ccRE_eRP.table.r.txt', 'r')
ccRE_eRP1 = {}

for records in ccRE_eRP0:
	tmp = [x.strip() for x in records.split('\t')]
	id = tmp[0]
	#if id in ccRE_eRP1:
	#	ccRE_eRP1[id].append(tmp)
	#else:
	ccRE_eRP1[id]=tmp

match = 0
ccRE_eRP_new = []
for records in ccRE1:
	if match%10000==0:
		print(match)
	if records in ccRE_eRP1:
		match = match+1
		ccRE_eRP_new.append(ccRE_eRP1[records])
	else:
		ccRE_eRP_new.append([records, '0', '0', '0', '0', '0', '0', '0', '0'])

### write 2d matrix
def write2d_array(array,output):
	r1=open(output,'w')
	for records in array:
		for i in range(0,len(records)-1):
			r1.write(str(records[i])+'\t')
		r1.write(str(records[len(records)-1])+'\n')
	r1.close()

ccRE_eRP_new = np.array(ccRE_eRP_new)
print(ccRE_eRP_new.shape)
write2d_array(ccRE_eRP_new, 'mouse_allChr_wTAD_022519.eRPs.txt')


