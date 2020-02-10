import os
import os.path
import numpy as np
from subprocess import call
from collections import Counter
import decimal

# create a new context for this task
ctx = decimal.Context()
ctx.prec = 20
################################################################################################
### float to str without scientific annotation
def float_to_str(f):
	"""
	Convert the given float to a string,
	without resorting to scientific notation
	"""
	d1 = ctx.create_decimal(repr(f))
	return format(d1, 'f')

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
allinfo = read2d_array('vision_rna.gene_ccRE.selected.all.reprod_count.WithNameSorted.txt', str)
#allinfo = read2d_array('test.txt', str)
if os.path.isfile('gene_bed_name.txt'):
	call('rm gene_bed_name.txt', shell=True)
if os.path.isfile('tss.allccREs.corabsmax.bed'):
	call('rm tss.allccREs.corabsmax.bed', shell=True)
if os.path.isfile('tss.allccREs.select.bed'):
	call('rm tss.allccREs.select.bed', shell=True)

current_gene = ''
for info in allinfo:
	gene = info[2]
	### get tss bed & tss ccREs
	if current_gene == '':
		print('start')
		current_gene = gene
		tss_chr = info[0]
		tss_start = info[1]
		start = float_to_str(int(float(info[1]))-1000000)
		end = float_to_str(int(float(info[1]))+1000000)
		if start >=0:
			tss2MB = np.array([[info[0], start, end]])
		else:
			tss2MB = np.array([[info[0], 0, end]])
		### tss bed
		write2d_array(tss2MB, 'tmp.tss.2MB.bed')
		### tss ccREs
		call('bedtools intersect -a vision_cres.sort.bed -b tmp.tss.2MB.bed -wa > tmp.tss.allccREs.bed', shell=True)
		#################
		#################
		### write 200bp ccREs
		ccRE_200bp = open('tmp.ccRE_200bp.txt','w')
		start_200bp = float_to_str(int(float(info[5])))
		end_200bp = float_to_str(int(float(info[5]))+200)
		if '1' in info[9]:
			ccRE_200bp.write(info[0] + '\t' + start_200bp + '\t' + end_200bp + '\t' + info[7] + '\t' + '1' + '\t' + str(info[9].count('1')) + '\n')
		else:
			ccRE_200bp.write(info[0] + '\t' + start_200bp + '\t' + end_200bp + '\t' + info[7] + '\t' + '0' + '\t' + '0' + '\n')
		#################
		#################
	elif gene == current_gene:
		### write 200bp ccREs
		start_200bp = float_to_str(int(float(info[5])))
		end_200bp = float_to_str(int(float(info[5]))+200)
		if '1' in info[9]:
			ccRE_200bp.write(info[0] + '\t' + start_200bp + '\t' + end_200bp + '\t' + info[7] + '\t' + '1' + '\t' + str(info[9].count('1')) + '\n')
		else:
			ccRE_200bp.write(info[0] + '\t' + start_200bp + '\t' + end_200bp + '\t' + info[7] + '\t' + '0' + '\t' + '0' + '\n')
		#################
		#################
	elif gene != current_gene:
		#################
		### close a previous gene
		#################
		ccRE_200bp.close()
		call('sort -k1,1 -k2,2n tmp.ccRE_200bp.txt > tmp.ccRE_200bp.sort.txt', shell=True)
		tmp_allccREs_linecount = len(open('tmp.tss.allccREs.bed').readlines())
		tmp_gene_bed_name = np.repeat([[tss_chr, tss_start, current_gene]], tmp_allccREs_linecount, axis=0)
		write2d_array(tmp_gene_bed_name, 'tmp.gene_bed_name.txt')
		call('bedtools map -a tmp.tss.allccREs.bed -b tmp.ccRE_200bp.sort.txt -c 4 -o absmax > tmp.tss.allccREs.corabsmax.bed', shell=True)
		call('bedtools map -a tmp.tss.allccREs.bed -b tmp.ccRE_200bp.sort.txt -c 6 -o max > tmp.tss.allccREs.select.bed', shell=True)
		### concatenate to all 
		call('cat tmp.gene_bed_name.txt >> gene_bed_name.txt', shell=True)
		call('cat tmp.tss.allccREs.corabsmax.bed >> tss.allccREs.corabsmax.bed', shell=True)
		call('cat tmp.tss.allccREs.select.bed >> tss.allccREs.select.bed', shell=True)
		#################
		### start a new gene
		#################
		print(gene)
		current_gene = gene
		tss_chr = info[0]	
		tss_start = info[1]
		start = float_to_str(int(float(info[1]))-1000000)
		end = float_to_str(int(float(info[1]))+1000000)
		if start >=0:
			tss2MB = np.array([[info[0], start, end]])
		else:
			tss2MB = np.array([[info[0], 0, end]])
		### tss bed
		write2d_array(tss2MB, 'tmp.tss.2MB.bed')
		### tss ccREs
		call('bedtools intersect -a vision_cres.bed -b tmp.tss.2MB.bed -wa > tmp.tss.allccREs.bed', shell=True)
		#################
		#################
		### write 200bp ccREs
		ccRE_200bp = open('tmp.ccRE_200bp.txt','w')
		start_200bp = float_to_str(int(float(info[5])))
		end_200bp = float_to_str(int(float(info[5]))+200)
		if '1' in info[9]:
			ccRE_200bp.write(info[0] + '\t' + start_200bp + '\t' + end_200bp + '\t' + info[7] + '\t' + '1' + '\t' + str(info[9].count('1')) + '\n')
		else:
			ccRE_200bp.write(info[0] + '\t' + start_200bp + '\t' + end_200bp + '\t' + info[7] + '\t' + '0' + '\t' + '0' + '\n')


#################
### get last info
#################
ccRE_200bp.close()
call('sort -k1,1 -k2,2n tmp.ccRE_200bp.txt > tmp.ccRE_200bp.sort.txt', shell=True)
tmp_allccREs_linecount = len(open('tmp.tss.allccREs.bed').readlines())
tmp_gene_bed_name = np.repeat([[tss_chr, tss_start, current_gene]], tmp_allccREs_linecount, axis=0)
write2d_array(tmp_gene_bed_name, 'tmp.gene_bed_name.txt')
call('bedtools map -a tmp.tss.allccREs.bed -b tmp.ccRE_200bp.sort.txt -c 4 -o absmax > tmp.tss.allccREs.corabsmax.bed', shell=True)
call('bedtools map -a tmp.tss.allccREs.bed -b tmp.ccRE_200bp.sort.txt -c 6 -o max > tmp.tss.allccREs.select.bed', shell=True)
### concatenate to all 
call('cat tmp.gene_bed_name.txt >> gene_bed_name.txt', shell=True)
call('cat tmp.tss.allccREs.corabsmax.bed >> tss.allccREs.corabsmax.bed', shell=True)
call('cat tmp.tss.allccREs.select.bed >> tss.allccREs.select.bed', shell=True)

### paste all
call('paste gene_bed_name.txt tss.allccREs.corabsmax.bed tss.allccREs.select.bed > all_converted.txt', shell=True)
call('cut -f1,2,3,5,6,7,11 all_converted.txt > all_converted.txt.tmp && mv all_converted.txt.tmp all_converted.txt', shell=True)







