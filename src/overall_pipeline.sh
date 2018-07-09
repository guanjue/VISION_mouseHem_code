###### read parameters from inputs
script_dir=$1
working_dir=$2
input_dir=$3
input_file_list=$4
overall_upper=$5
overall_lower=$6
select_method=$7
user_given_global_ref=$8


###### covert reads count to NB p-value
while read LINE
do
	sig1=$(echo "$LINE" | awk '{print $1}')
	sig2=$(echo "$LINE" | awk '{print $2}')
	echo $sig1 $sig2
	if time Rscript $script_dir'negative_binomial_p_2r_bgadj_bayes.R' $sig1 $input_dir $sig2 $input_dir $sig1; then echo 'covert reads count to NB p-value DONE'; else echo 'ERROR: covert reads count to NB p-value' && exit 1; fi
done < $input_file_list



###### convert nbp to Fisher's method merged p-value
### extrac cell type mark list
ls *.nbp_2r_bgadj.txt | awk -F '.' -v OFS='\t' '{print $1"."$2}' | sort -u > cell_marker_list.txt
ls *.nbp_2r_bgadj.txt | awk -F '.' -v OFS='\t' '{print $2}' | sort -u > mark_list.txt
ls *.nbp_2r_bgadj.txt | awk -F '.' -v OFS='\t' '{print $1}' | sort -u > cell_list.txt

### move data NB p-value data into nbp folder
if [ -d $working_dir'nbp/' ]; then echo $working_dir'nbp/' exist; else mkdir $working_dir'nbp/'; fi
mv *.nbp_2r_bgadj.txt $working_dir'nbp/'
mv *.mvsp.txt $working_dir'nbp/'
### get Fisher's method merged pval
for cm in $(cat cell_marker_list.txt)
do
	echo $cm
	if time Rscript $script_dir'fisher_pval.R' $cm '.nbp_2r_bgadj.txt' $working_dir'nbp/' 100; then echo 'get Fisher method merged pval DONE'; else echo 'ERROR: get Fishers method merged pval' && exit 1; fi
done



###### select reference dataset for pknorm
for mk in $(cat mark_list.txt)
do
	echo $mk
	ls *$mk*.frip_snr.txt > $mk'.file_list.txt'
	if time Rscript $script_dir'get_mk_ref.R' $mk'.file_list.txt' $select_method $mk'.ref_frip.txt'; then echo 'select reference dataset for pknorm'; else echo 'ERROR: select reference dataset for pknorm' && exit 1; fi
done
### select top reference dataset for cross mark pknorm
if time Rscript $script_dir'get_top_ref.R' '.ref_frip.txt' $select_method $working_dir cross_mark_ref_list.txt $user_given_global_ref; then echo 'select top reference dataset for cross mark pknorm DONE'; else echo 'ERROR: select top reference dataset for cross mark pknorm' && exit 1; fi



###### pknorm normalize reference datasets of all marks
while read LINE
do
	sig1=$(echo "$LINE" | awk '{print $1}')
	sig2=$(echo "$LINE" | awk '{print $2}')
	sig2_celltype=$(echo "$LINE" | awk '{print $2}' | awk -F '.' -v OFS='\t' '{print $1"_"$2}')
	upperlim=100
	lowerlim=0
	echo $sig1 
	echo $sig2
	echo $sig2_celltype
	### set upper limit
	cat $sig1 | awk -F '\t' -v OFS='\t' -v ul=$upperlim '{if ($1>=ul) print ul; else print $1}' > $sig1'.upperlim.txt'
	cat $sig2 | awk -F '\t' -v OFS='\t' -v ul=$upperlim '{if ($1>=ul) print ul; else print $1}' > $sig2'.upperlim.txt' 
	### peak norm
	if time python $script_dir'peaknorm_rotate_log_ref_mean.py' -n 500000 -a $sig1'.upperlim.txt' -b $sig2'.upperlim.txt' -u $upperlim -l $lowerlim; then echo 'pknorm normalize reference DONE'; else echo 'ERROR: pknorm normalize reference' && exit 1; fi
	### rm tmp files
	rm $sig1'.upperlim.txt'
	rm $sig2'.upperlim.txt'
done < cross_mark_ref_list.txt.info.txt
### move ref norm files into ref_info folder
if [ -d $working_dir'ref_info/' ]; then echo $working_dir'ref_info/' exist; else mkdir $working_dir'ref_info/'; fi
mv *.pknorm.scatterplot.png $working_dir'ref_info/'
mv *.scatterplot.png $working_dir'ref_info/'
mv *.ref.info.txt $working_dir'ref_info/'



###### pknorm across datasets with the same mark
for mk in $(cat mark_list.txt)
do
	echo $mk
	while read LINE
	do
		sig1=$(echo "$LINE" | awk '{print $1}')
		sig2=$(echo "$LINE" | awk '{print $2}')
		sig2_celltype=$(echo "$LINE" | awk '{print $2}' | awk -F '.' -v OFS='\t' '{print $1"_"$2}')
		upperlim=100
		lowerlim=0
		echo $sig1 
		echo $sig2
		echo $sig2_celltype
		### set upper limit
		cat $sig1 | awk -F '\t' -v OFS='\t' -v ul=$upperlim '{if ($1>=ul) print ul; else print $1}' > $sig1'.upperlim.txt'
		cat $sig2 | awk -F '\t' -v OFS='\t' -v ul=$upperlim '{if ($1>=ul) print ul; else print $1}' > $sig2'.upperlim.txt' 
		### peak norm
		if time python $script_dir'peaknorm_rotate_log_z_mean.py' -n 500000 -a $sig1'.upperlim.txt' -b $sig2'.upperlim.txt' -u $upperlim -l $lowerlim; then echo 'pknorm across datasets DONE'; else echo 'ERROR: pknorm across datasets' && exit 1; fi
		### rm tmp files
		rm $sig1'.upperlim.txt'
		rm $sig2'.upperlim.txt'
	done < $mk'.ref_frip.txt.info.txt'
done
### move pknorm files into pknorm_sig folder
if [ -d $working_dir'pknorm_info/' ]; then echo $working_dir'pknorm_info/' exist; else mkdir $working_dir'pknorm_info/'; fi
mv *.pknorm.scatterplot.png $working_dir'pknorm_info/'
mv *.scatterplot.png $working_dir'pknorm_info/'
mv *.info.txt $working_dir'pknorm_info/'



###### mv PKnorm normalized signal files & unnormalized signal files into *_sig folders
### pknorm signal
if [ -d $working_dir'pknorm_sig/' ]; then echo $working_dir'pknorm_sig/' exist; else mkdir $working_dir'pknorm_sig/'; fi
mv *.pknorm.txt $working_dir'pknorm_sig/'
### ref pknorm signal
if [ -d $working_dir'pknorm_ref_sig/' ]; then echo $working_dir'pknorm_ref_sig/' exist; else mkdir $working_dir'pknorm_ref_sig/'; fi
mv *.pknorm.ref.txt $working_dir'pknorm_ref_sig/'
### ref frip & snp
mv *.ref_frip.txt $working_dir'ref_info/'
### fisher pvalue signal without normalization
if [ -d $working_dir'fisherp/' ]; then echo $working_dir'fisherp/' exist; else mkdir $working_dir'fisherp/'; fi
mv *.fisher_p.txt $working_dir'fisherp/'
mv *.frip_snr.txt $working_dir'fisherp/'



###### set limit for signals
for filename in $(cat cell_marker_list.txt)
do
	echo $filename
	if cat $working_dir'pknorm_sig/'$filename'.pknorm.txt' | awk -F '\t' -v OFS='\t' -v ul=$overall_upper -v ll=$overall_lower '{if ($1<ll) print ll; else if ($1>ul) print ul; else print $1}' > $filename'.pknorm.'$overall_lower'_'$overall_upper'.txt'; then echo 'set limit for signals DONE'; else echo 'ERROR: set limit for signals' && exit 1; fi
done
### mv to the output folder
if [ -d $working_dir'pknorm_'$overall_lower'_'$overall_upper'_sig/' ]; then echo $working_dir'pknorm_'$overall_lower'_'$overall_upper'_sig/' exist; else mkdir $working_dir'pknorm_'$overall_lower'_'$overall_upper'_sig/'; fi
mv *'.pknorm.'$overall_lower'_'$overall_upper'.txt' $working_dir'pknorm_'$overall_lower'_'$overall_upper'_sig/'



### list files
if [ -d $working_dir'list_files/' ]; then echo $working_dir'list_files/' exist; else mkdir $working_dir'list_files/'; fi
mv *list.txt $working_dir'list_files/'





