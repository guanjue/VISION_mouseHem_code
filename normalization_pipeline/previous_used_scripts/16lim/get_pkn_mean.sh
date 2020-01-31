#PBS -l nodes=4:ppn=8
#PBS -l walltime=40:00:00
#PBS -j oe
#PBS -A yzz2_c_t_sc_default

module load gcc
module load python/2.7.8

cd /storage/home/gzx103/scratch/vision/5end/pknorm_16lim/pcor_100lim_mean/

#sleep 600
rm *.pknorm.txt
rm *.bw
rm *.pkn.txt
rm *.pkn16.txt
for mark in $(cat mark_list.txt)
do
	echo $mark
	cp '/storage/home/gzx103/scratch/vision/5end/pknorm_16lim/'$mark'_r_log_mm/'*.pknorm.txt ./
	ls *.pknorm.txt > pknorm_list.tmp.txt
	for file in $(cat pknorm_list.tmp.txt)
	do
		ct=$(echo "$file" | awk -F '.' '{print $1}')
		echo $ct
		mv $ct'.pknorm.txt' $ct'.'$mark'.pkn.txt'
		cat $ct'.'$mark'.pkn.txt' | awk -F '\t' -v OFS='\t' '{if ($1>=16) print 16; else print $1}' > $ct'.'$mark'.pkn16.txt'
	done
done

mkdir pknorm_16lim_ref1mo_sample1mo_correct
cp *.pkn16.txt pknorm_16lim_ref1mo_sample1mo_correct
cp -r pknorm_16lim_ref1mo_sample1mo_correct /storage/home/gzx103/group/projects/vision/
cd /storage/home/gzx103/group/projects/vision/
chmod 777 pknorm_16lim_ref1mo_sample1mo_correct
cd -


ls pknorm_16lim_ref1mo_sample1mo/*.pkn16.txt > pkn16.list.txt
### convert sig to bigwig
for filename in $(cat pkn16.list.txt)
do
        echo $filename
        ct_mark=$(echo "$filename" | awk -F '/' '{print $2}'| awk -F '.' '{print $1"."$2}')
        echo $ct_mark
        ### get bedgraph
        paste ~/group/projects/vision/merged_input/200_noblack.11_22_2017.bed 'pknorm_16lim_ref1mo_sample1mo/'$ct_mark'.pkn16.txt' | awk -F '\t' -v OFS='\t' '{print $1,$2,$3,$4}' | sort -k1,1 -k2,2n> $ct_mark'.1mo.bedgraph'
        ### bedGraphToBigWig
        /storage/home/gzx103/group/projects/vision/input_norm/bedGraphToBigWig $ct_mark'.1mo.bedgraph' /storage/home/gzx103/group/projects/vision/input_norm/mm10.chrom.sizes $ct_mark'.1mo.16lim.bw'
        ### rm bedgraph
        rm $ct_mark'.1mo.bedgraph'
done

rm -r /storage/home/gzx103/scratch/vision/5end/pknorm_16lim/pcor_100lim_mean/bw_pkn16
mkdir /storage/home/gzx103/scratch/vision/5end/pknorm_16lim/pcor_100lim_mean/bw_pkn16
mv *.1mo.16lim.bw /storage/home/gzx103/scratch/vision/5end/pknorm_16lim/pcor_100lim/bw_pkn16/

mkdir pknorm_16lim_ref1mo_sample1mo_bw
cp *.1mo.bw pknorm_16lim_ref1mo_sample1mo_bw/


ls *.pkn.txt > pkn_list.txt
time Rscript /storage/home/gzx103/scratch/vision/5end/pknorm_16lim/pcor_100lim/cor_heatmap.R

mkdir pcor_pair
cp pk_cor.matrix.txt pcor_pair

mkdir pcor_pair_sp
cp pk_cor.matrix.sp.txt pcor_pair_sp

cd pcor_pair


