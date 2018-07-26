import numpy as np
from subprocess import call

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

def plot_ven(file_list_file):
	file_list = read2d_array(file_list_file, 'str')

	output = open(file_list_file+'.venn.R', 'w')

	output.write('library(venneuler)\n')
	output.write('vc = venneuler(c(')
	for i in range(0, len(file_list)):
		pk = read2d_array(file_list[i][0], 'str')
		print(file_list[i][0])
		pk_num = pk.shape[0]
		print(pk_num)
		label_all = file_list[i][0].split('.')[0].split('_')
		label = ''
		for l in label_all:
			if (l!='0') and (l!='X'):
				if label=='':
					label = label+l
				else:
					label = label+'&'+l
		print(label)
		if i != len(file_list)-1:
			output.write('\'' + label + '\'' + '=' + str(pk_num) + ', ')
		else:
			output.write('\'' + label + '\'' + '=' + str(pk_num))


	output.write('))\n')
	output.write('pdf(\''+file_list_file+'.pdf\')\n')
	output.write('plot(vc)\n')
	output.write('dev.off()')
	output.close()

	call('Rscript ' +file_list_file+'.venn.R', shell=True)

############################################################################

import getopt
import sys
def main(argv):
	try:
		opts, args = getopt.getopt(argv,"hl:")
	except getopt.GetoptError:
		print 'time python plot_ven.py -l VAF_list.txt'
		sys.exit(2)

	for opt,arg in opts:
		if opt=="-h":
			print 'time python plot_ven.py -l VAF_list.txt'		
		elif opt=="-l":
			file_list=str(arg.strip())

	plot_ven(file_list)

if __name__=="__main__":
	main(sys.argv[1:])


