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

################################################################################################
### get bedgraph
def get_bedgraph(input, output):
	### read input
	allinfo = read2d_array(input, str)
	### check selection
	selinfo = allinfo[allinfo[:,8]=='1',:]
	### count selection
	selc = []
	for sel in selinfo[:,9]:
		selc.append([sel.count('1')])
	### reshape
	selc = np.array(selc)
	selc.reshape(selc.shape[0],1)
	print(selc.shape)
	### concatenate
	chrom = selinfo[:,0].reshape(selc.shape[0],1)
	start = selinfo[:,5].reshape(selc.shape[0],1)
	end = ((selinfo[:,5].astype(np.int))+200).reshape(selc.shape[0],1)
	cor = selinfo[:,7].reshape(selc.shape[0],1)
	print(chrom.shape)
	print(start.shape)
	print(end.shape)
	print(cor.shape)
	output_mat1 = np.concatenate((chrom, start, end, cor), 1)
	output_mat2 = np.concatenate((chrom, start, end, selc), 1)
	### write output
	write2d_array(output_mat1, output+'.cor.bedgraph')
	write2d_array(output_mat2, output+'.repr.bedgraph')



###########################################
# time python extract_bedgraph.py -i Hba-a1.ENSMUSG00000069919.7.allinfo.txt -o Hba-a1.ENSMUSG00000069919.7
import getopt
import sys

def main(argv):
	try:
		opts, args = getopt.getopt(argv,"hi:o:")
	except getopt.GetoptError:
		print 'python extract_statepair_position.py -i sample1.tab -j sample2.tab -o output_file.txt'
		sys.exit(2)

	for opt,arg in opts:
		if opt=="-h":
			print 'python extract_statepair_position.py -i sample1.tab -j sample2.tab -o output_file.txt'
			sys.exit()
		elif opt=="-i":
			input=str(arg.strip())
		elif opt=="-o":
			output=str(arg.strip())
	get_bedgraph(input, output)

if __name__=="__main__":
	main(sys.argv[1:])
