import numpy as np

def ct_dist_ideas_state_mean_signal(ideas_state_matrix, ideas_meansig_matrix, output_name):
	################################################################################################
	### write 2d matrix
	def write2d_array(array,output):
		r1=open(output,'w')
		for records in array:
			for i in range(0,len(records)-1):
				r1.write(str(records[i])+'\t')
			r1.write(str(records[len(records)-1])+'\n')
		r1.close()
	### random seed
	np.random.seed(2018)
	############
	###### read whole genome IDEAS state in each cell type
	ideas_state_data0 = open(ideas_state_matrix, 'r')
	ct_name = ideas_state_data0.readline().split(' ')[4:-1]
	print(ct_name)
	###### get IDEAS state matrix
	ideas_state_data = []
	readline_i = 0
	for bin_state_vector in ideas_state_data0:
		readline_i = readline_i+1
		if readline_i%10000 == 0:
			print(readline_i)
		bin_state_vector_tmp = bin_state_vector.split(' ')[4:-1]
		ideas_state_data.append(bin_state_vector_tmp)

	ideas_state_data = np.array(ideas_state_data)
	idx = np.random.randint(ideas_state_data.shape[0], size=100000) 
	ideas_state_data = ideas_state_data[idx,:]
	ideas_state_data0.close()

	print(np.sum(ideas_state_data[:,0]!=ideas_state_data[:,1]))
	############
	###### read IDEAS state mean signal
	ideas_state_meansig0 = open(ideas_meansig_matrix, 'r')
	mark_list = ideas_state_meansig0.readline().split()
	###### generate IDEAS state mean signal dict
	ideas_state_meansig_dict = {}
	for state_meansig_vector in ideas_state_meansig0:
		state_meansig_vector_tmp = state_meansig_vector.split()
		state_tmp = state_meansig_vector_tmp[0]
		signal_tmp = np.array(state_meansig_vector_tmp[1:], dtype=float)
		ideas_state_meansig_dict[state_tmp] = signal_tmp

	ideas_state_meansig0.close()


	ct_num = ideas_state_data.shape[1]
	Euclidean_dist_matrix = np.zeros((ct_num, ct_num))
	Euclidean_dist_matrix_dict = {}
	for mark_i in range(0,8):
		Euclidean_dist_matrix_dict[mark_i] = np.zeros((ct_num, ct_num))

	for i in range(0,ct_num):
		for j in range(0,ct_num):
			if i > j:
				print(str(i)+'_'+str(j))
				ct1_ideas = ideas_state_data[:,i]
				ct2_ideas = ideas_state_data[:,j]
				ct1_ideas_signal_matrix = []
				ct2_ideas_signal_matrix = []
				for cts1, cts2 in zip(ct1_ideas, ct2_ideas):
					ct1_ideas_signal_matrix.append(ideas_state_meansig_dict[cts1])
					ct2_ideas_signal_matrix.append(ideas_state_meansig_dict[cts2])
				############
				###### calculate mv Euclidean distance
				Euclidean_dist_square = 0.0
				Euclidean_dist_square_dict = {}
				for mark_i in range(0,8):
					Euclidean_dist_square_dict[mark_i] = 0.0

				for sig1, sig2 in zip(ct1_ideas_signal_matrix, ct2_ideas_signal_matrix):
					Euclidean_dist_square = Euclidean_dist_square + np.sum(np.square(sig1-sig2))
					for mark_i in range(0,8):
						Euclidean_dist_square_dict[mark_i] = Euclidean_dist_square_dict[mark_i] + (np.square(sig1-sig2)[mark_i])

				Euclidean_dist = np.sqrt(Euclidean_dist_square)
				print(Euclidean_dist)
				Euclidean_dist_matrix[i,j] = Euclidean_dist
				Euclidean_dist_matrix[j,i] = Euclidean_dist
				for mark_i in range(0,8):
					Euclidean_dist_mark_i = np.sqrt(Euclidean_dist_square_dict[mark_i])
					print(Euclidean_dist_mark_i)
					Euclidean_dist_matrix_dict[mark_i][i,j] = Euclidean_dist_mark_i
					Euclidean_dist_matrix_dict[mark_i][j,i] = Euclidean_dist_mark_i


	Euclidean_dist_matrix = np.array(Euclidean_dist_matrix)
	write2d_array(Euclidean_dist_matrix, output_name)

	for mark_i in range(0,8):
		Euclidean_dist_matrix_mark_i = np.array(Euclidean_dist_matrix_dict[mark_i])
		write2d_array(Euclidean_dist_matrix_mark_i, mark_list[i] + '_' + output_name)






############################################################################
#time python ct_dist_ideas_state_mean_signal.py -i pknorm_2_16lim_ref1mo_0424_lesshet.state -s ideas_state_mean_signal.txt -o IDEAS_state_meansignal_Euclidean_dist_matrix.txt

import getopt
import sys
def main(argv):
	try:
		opts, args = getopt.getopt(argv,"hi:s:o:")
	except getopt.GetoptError:
		print 'time python ct_dist_ideas_state_mean_signal.py -i ideas_state_matrix -s ideas_meansig_matrix -o output_name'
		sys.exit(2)

	for opt,arg in opts:
		if opt=="-h":
			print 'time python ct_dist_ideas_state_mean_signal.py -i ideas_state_matrix -s ideas_meansig_matrix -o output_name'
			sys.exit()
		elif opt=="-i":
			ideas_state_matrix=str(arg.strip())			
		elif opt=="-s":
			ideas_meansig_matrix=str(arg.strip())					
		elif opt=="-o":
			output_name=str(arg.strip())		


	ct_dist_ideas_state_mean_signal(ideas_state_matrix, ideas_meansig_matrix, output_name)

if __name__=="__main__":
	main(sys.argv[1:])


