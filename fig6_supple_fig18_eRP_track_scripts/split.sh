ct_list_file='/storage/home/gzx103/group/projects/vision/pknorm_2_16lim_ref1mo_sample1mo_0424/get_eRP_tracks/ct_list.txt'
output_folder='/storage/home/gzx103/group/projects/vision/pknorm_2_16lim_ref1mo_sample1mo_0424/get_eRP_tracks/'
get_8tracks_file='/storage/home/gzx103/group/projects/vision/pknorm_2_16lim_ref1mo_sample1mo_0424/get_eRP_tracks/get_8tracks.R'

for ct in $(cat $ct_list_file)
do
	### split for speed up
	mkdir $ct'_folder'
	for i in {0..61}
	do
		echo $i
		rownum=$((i * 10000 + 1))
		tail -n+$rownum 'ccRE.'$ct'.state.2.pksort.bed' | head -10000 > 'ccRE.'$ct'.state.2.'$i'.bed'
		mv 'ccRE.'$ct'.state.2.'$i'.bed' $ct'_folder'
	done
	### get IDEAS state percentage tracks
	for i in {0..61}
	do
		echo \#PBS -l nodes=1:ppn=8 > 'PBS.'$ct$i'.sh'
		echo \#PBS -l walltime=10:00:00 >> 'PBS.'$ct$i'.sh'
		echo \#PBS -j oe >> 'PBS.'$ct$i'.sh'
		echo \#PBS -A yzz2_e_g_sc_default >> 'PBS.'$ct$i'.sh'
		echo \#PBS -l pmem=16gb >> 'PBS.'$ct$i'.sh'
		echo module load gcc/5.3.1 >> 'PBS.'$ct$i'.sh'
		echo module load python/2.7.14-anaconda5.0.1 >> 'PBS.'$ct$i'.sh'
		echo module load bedtools >> 'PBS.'$ct$i'.sh'
		echo cd $output_folder$ct'_folder' >> 'PBS.'$ct$i'.sh' 
		echo time Rscript $get_8tracks_file 'ccRE.'$ct'.state.2.'$i'.bed' 'ccRE.'$ct'.state.2.'$i'.statepercent.bed' >> 'PBS.'$ct$i'.sh'
		qsub 'PBS.'$ct$i'.sh'
	done
done
