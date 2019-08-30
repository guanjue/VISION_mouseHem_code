'''
##################
### download bigbed files from codex
for file in $(cat file_list.txt)
do
	echo $file
	wget 'http://codex.stemcells.cam.ac.uk/data/bb/mm10/'$file
done
mkdir bb_files
mv *bb bb_files

##################
### convert bigbed to bed
mkdir bed_files_mm10
for file in $(cat file_list.txt)
do
	echo $file
	time ~/group/software/ucsc/bigBedToBed bb_files/$file bed_files_mm10/$file'.bed'
	time sort -k1,1 -k2,2n bed_files_mm10/$file'.bed' > bed_files_mm10/$file'.sort.bed'
done
'''

##################
### jarkard index
for file in $(cat ../file_list.sub.txt)
do
	echo $file
	i=$(echo $file | awk -F '.' '{print $1}')
	echo $i
	mkdir 'jaccard_index_'$i
	time sort -k1,1 -k2,2n '/storage/home/gzx103/group/projects/vision/snapshot18_reproduce_0_16lim_all_NB88/index_set_bed/'$file > IS_tmp.sort.bed
	### intersect
	for peak_file in $(cat ../file_list.sub.txt)
	do
		echo $peak_file
		time bedtools jaccard -a IS_tmp.sort.bed -b ../bed_files_mm10/$peak_file'.sort.bed' > 'jaccard_index_'$i'/IS'$i'.'$peak_file'.JI.txt' 
	done
done


##################
### Jaccard index matrix
mkdir JI_matrix_folder 
for file in $(cat ../file_list.sub.txt)
do
	echo $file
	i=$(echo $file | awk -F '.' '{print $1}')
	echo $i
	for peak_file in $(cat ../file_list.sub.txt)
	do
		echo $peak_file
		time tail -n+2 'jaccard_index_'$i'/IS'$i'.'$peak_file'.JI.txt' >> 'JI_matrix_'$i'.txt'
	done
	mv 'JI_matrix_'$i'.txt' JI_matrix_folder
done






