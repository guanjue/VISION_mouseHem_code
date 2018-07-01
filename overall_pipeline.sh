

script_dir=/storage/home/gzx103/group/software/signorm/bin/

script_dir=/storage/home/gzx103/group/software/vision_mouse/src/
working_dir=/storage/home/gzx103/scratch/vision/5end/


###### covert reads count to NB p-value
while read LINE
do
	sig1=$(echo "$LINE" | awk '{print $1}')
	sig2=$(echo "$LINE" | awk '{print $2}')
	echo $sig1 $sig2
	time Rscript $script_dir'negative_binomial_p_2r_bgadj.R' $sig1 $sig2 $sig1'.nbp_bgadj.txt'
done < info_table_all.rc2nbp.txt


###### convert nbp to Fisher's method merged p-value
### extrac cell type marker list
ls *.nbp_bgadj.txt.nbp_2r_bgadj.txt | awk -F '.' -v OFS='\t' '{print $1"."$2}' > cell_marker_list.txt
if [ -d $working_dir'nbp/' ]; then rm -Rf $WORKING_DIR; else mkdir $working_dir'nbp/'; fi
mv *.nbp_bgadj.txt.nbp_2r_bgadj.txt $working_dir'nbp/'
### get Fisher's method combined pval
for cm in $(cat cell_marker_list.txt)
do
	echo $cm
	Rscript $script_dir'fisher_pval.R' $cm '.nbp_bgadj.txt.nbp_2r_bgadj.txt' $working_dir'nbp/' 100
done


###### select reference dataset for pknorm
for mk in $(cat mark_list.txt)
do
	echo $mk
	ls *$mk*.mvsp.txt > $mk'.file_list.txt'
	time Rscript $script_dir'get_mk_ref.R' $mk'.file_list.txt' $mk'.ref_frip.txt'
done


###### pknorm normalize reference datasets of all marks
while read LINE
do
	sig1=$(echo "$LINE" | awk '{print $1}')
	sig1_col=$(echo "$LINE" | awk '{print $2}')
	sig2=$(echo "$LINE" | awk '{print $3}')
	sig2_col=$(echo "$LINE" | awk '{print $4}')
	sig2_celltype=$(echo "$LINE" | awk '{print $3}' | awk -F '.' -v OFS='\t' '{print $1"_"$2}')
	upperlim=100
	lowerlim=0
	echo $sig1 
	echo $sig1_col 
	echo $sig2
	echo $sig2_col
	echo $sig2_celltype
	### set upper limit
	cat $sig1 | awk -F '\t' -v OFS='\t' -v ul=$upperlim '{if ($1>=ul) print ul; else print $1}' > $sig1'.upperlim.txt'
	cat $sig2 | awk -F '\t' -v OFS='\t' -v ul=$upperlim '{if ($1>=ul) print ul; else print $1}' > $sig2'.upperlim.txt' 
	### peak norm
	time python $script_dir'peaknorm_rotate_log_ref_mean.py' -w /storage/home/gzx103/scratch/vision/5end/pknorm_16lim/200_noblack.11_22_2017.bed -p $mark'TablePeaksFiltered.txt' -n 500000 -a $sig1_col -b $sig1'.upperlim.txt' -c $sig2_col -d $sig2'.upperlim.txt' -u $upperlim -l $lowerlim
done < raw_sig_list.txt







###### pknorm across datasets with the same mark




###### set limit for signals

