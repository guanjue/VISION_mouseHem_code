

script_dir=/storage/home/gzx103/group/software/signorm/bin/
working_dir=/storage/home/gzx103/scratch/vision/5end/


###### covert reads count to NB p-value
while read LINE
do
	sig1=$(echo "$LINE" | awk '{print $1}')
	sig2=$(echo "$LINE" | awk '{print $2}')
	echo $sig1 $sig2
	time Rscript $script_folder'negative_binomial_p_2r_bgadj.R' $sig1 $sig2 $sig1'.nbp_bgadj.txt'
done < info_table_all.rc2nbp.txt


###### convert nbp to Fisher's method merged p-value
### extrac cell type marker list
ls *.nbp_bgadj.txt.nbp_2r_bgadj.txt | awk -F '.' -v OFS='\t' '{print $1"."$2}' > cell_marker_list.txt
mkdir $working_dir'nbp/'
mv *.nbp_bgadj.txt.nbp_2r_bgadj.txt $working_dir'nbp/'
### get Fisher's method combined pval
for cm in $(cat cell_marker_list.txt)
do
	echo $cm
	Rscript $script_folder'fisher_pval.R' $cm '.nbp_bgadj.txt.nbp_2r_bgadj.txt' $working_dir'nbp/' 100
done


