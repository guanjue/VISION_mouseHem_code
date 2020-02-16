#PBS -l nodes=1:ppn=4
#PBS -l walltime=20:00:00
#PBS -j oe
#PBS -A open

module load gcc
module load python/2.7.8

##################################
script_folder='/storage/home/gzx103/group/software/signorm/bin/'
analysis_folder='/storage/home/gzx103/group/projects/vision/merged_input/'

input_folder='/storage/home/gzx103/group/projects/vision/merged_input/'
output_folder_t_r_file='/storage/home/gzx103/group/projects/vision/merged_input/'
output_folder_normed_sig_file='/storage/home/gzx103/group/projects/vision/merged_input/'
##################################
cd $analysis_folder
##################################
##################################
#rm raw_input_read.txt
#while read LINE
#do
#        sig1=$(echo "$LINE" | awk '{print $1}')
#        sig2=$(echo "$LINE" | awk '{print $2}')
#	paste -sd+ $sig1 | bc >> raw_input_read.txt
#done < $analysis_folder'info_table_input.txt'
##################################
### run pipeline
while read LINE
do
	sig1=$(echo "$LINE" | awk '{print $1}')
	sig2=$(echo "$LINE" | awk '{print $2}')
	echo $sig1 $sig2
	### get ncis t_r matrix
	time python $script_folder'get_ncis_t_a_b.py' -i $input_folder$sig1 -j $input_folder$sig2 -o $output_folder_t_r_file$sig1'_vs_'$sig2'.txt'
	### get scale factor and normalize x-axis signal
	time Rscript $script_folder'signorm.R' $input_folder$sig1 $input_folder$sig2 $output_folder_t_r_file$sig1'_vs_'$sig2'.txt' $output_folder_normed_sig_file$sig1'.totalmean.txt' meanvar BinSeg 5 polynorm 50 1000000 1 2017 0 1000000 1 1 $script_folder 1 linear
done < $analysis_folder'info_table_input.txt'
##################################
### merge normed input
paste $output_folder_normed_sig_file*'.totalmean.txt' | awk -F '\t' -v OFS='\t' '{sum=0; for(i=1; i<=NF; i++){sum+=$i}; sum/=NF; print sum}' > $output_folder_normed_sig_file'merged_normed_input_mean.txt'
paste $output_folder_normed_sig_file*'.totalmean.txt' | awk -F '\t' -v OFS='\t' '{sum=0; for(i=1; i<=NF; i++){sum+=$i}; print sum}' > $output_folder_normed_sig_file'merged_normed_input_sum.txt'

### rounding $output_folder_normed_sig_file'merged_normed_input.txt
cat $output_folder_normed_sig_file'merged_normed_input_mean.txt' | awk -F '\t' -v OFS='\t' '{printf "%1.0f\n", $1}' > $output_folder_normed_sig_file'merged_normed_input_mean.rounding.txt'
cp $output_folder_normed_sig_file'merged_normed_input_mean.rounding.txt' $output_folder_normed_sig_file'merged_normed_input.rounding.txt'
cat $output_folder_normed_sig_file'merged_normed_input_sum.txt' | awk -F '\t' -v OFS='\t' '{printf "%1.0f\n", $1}' > $output_folder_normed_sig_file'merged_normed_input_sum.rounding.txt'

paste ~/group/projects/vision/merged_input/200_noblack.11_22_2017.bed merged_normed_input.rounding.txt | awk -F '\t' -v OFS='\t' '{print $1,$2,$3,$4}' | sort -k1,1 -k2,2n> merged_normed_input.rounding.bedgraph

./bedGraphToBigWig merged_normed_input.rounding.bedgraph /storage/home/gzx103/group/projects/vision/input_norm/mm10.chrom.sizes merged_normed_input.rounding.bw


paste ~/group/projects/vision/merged_input/200_noblack.11_22_2017.bed 200_noblack.5bins_bgsig_mean.round.txt | awk -F '\t' -v OFS='\t' '{print $1,$2,$3,$4}' | sort -k1,1 -k2,2n> 200_noblack.5bins_bgsig_mean.round.bedgraph

/storage/home/gzx103/group/projects/vision/input_norm/bedGraphToBigWig 200_noblack.5bins_bgsig_mean.round.bedgraph /storage/home/gzx103/group/projects/vision/input_norm/mm10.chrom.sizes 200_noblack.5bins_bgsig_mean.round.bw

rm 200_noblack.5bins_bgsig_mean.round.bedgraph


paste ~/group/projects/vision/merged_input/200_noblack.11_22_2017.bed 200_noblack.5bins_bgsig_sum.round.txt | awk -F '\t' -v OFS='\t' '{print $1,$2,$3,$4}' | sort -k1,1 -k2,2n> 200_noblack.5bins_bgsig_sum.round.bedgraph

/storage/home/gzx103/group/projects/vision/input_norm/bedGraphToBigWig 200_noblack.5bins_bgsig_sum.round.bedgraph /storage/home/gzx103/group/projects/vision/input_norm/mm10.chrom.sizes 200_noblack.5bins_bgsig_sum.round.bw

rm 200_noblack.5bins_bgsig_sum.round.bedgraph

