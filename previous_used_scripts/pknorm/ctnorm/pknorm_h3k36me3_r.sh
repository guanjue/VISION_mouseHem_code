
#PBS -l nodes=2:ppn=8
#PBS -l walltime=10:00:00
#PBS -j oe
#PBS -A yzz2_c_t_sc_default

module load gcc
module load python/2.7.8


mark=h3k36me3
cd '/storage/home/gzx103/scratch/vision/5end/pknorm_16lim/'$mark'_r_log_mm'

cp '/storage/home/gzx103/scratch/vision/5end/fisher_p_100lim/'*'.'$mark'rep.fisher_p.txt' ./

while read LINE
do
	sig1=$(echo "$LINE" | awk '{print $1}')
	sig1_col=$(echo "$LINE" | awk '{print $2}')
	sig2=$(echo "$LINE" | awk '{print $3}')
	sig2_col=$(echo "$LINE" | awk '{print $4}')
	sig2_celltype=$(echo "$LINE" | awk '{print $3}' | awk -F '.' -v OFS='\t' '{print $1}')
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
	time python /storage/home/gzx103/group/software/signorm/bin/peaknorm_rotate_log_z_mean.py -w /storage/home/gzx103/scratch/vision/5end/pknorm_16lim/200_noblack.11_22_2017.bed -p $mark'TablePeaksFiltered.txt' -n 500000 -a $sig1_col -b $sig1'.upperlim.txt' -c $sig2_col -d $sig2'.upperlim.txt' -u $upperlim -l $lowerlim
	### get bed files
	paste /storage/home/gzx103/scratch/vision/5end/pknorm_16lim/200_noblack.11_22_2017.bed $sig2_celltype'.'$mark'rep.fisher_p.txt' | awk -F '\t' -v OFS='\t' '{print $1, $2, $3, ".", $4}' > $sig2_celltype'.'$mark'rep.fisher_p.bed'
	paste /storage/home/gzx103/scratch/vision/5end/pknorm_16lim/200_noblack.11_22_2017.bed $sig2_celltype'.pknorm.txt' | awk -F '\t' -v OFS='\t' '{print $1, $2, $3, ".", $4}' > $sig2_celltype'.pknorm.bed'
	cat $sig2_celltype'.wg.bed' | awk -F '\t' -v OFS='\t' '{if ($4!=0) print $1, $2, $3, ".", 1}' > $sig2_celltype'.pk.bed'

	paste ~/group/projects/vision/merged_input/200_noblack.11_22_2017.bed $sig2_celltype'.pknorm.txt' | awk -F '\t' -v OFS='\t' '{print $1,$2,$3,$4}' | sort -k1,1 -k2,2n> $sig2_celltype'.bedgraph'
	### bedGraphToBigWig
	/storage/home/gzx103/group/projects/vision/input_norm/bedGraphToBigWig $sig2_celltype'.bedgraph' /storage/home/gzx103/group/projects/vision/input_norm/mm10.chrom.sizes $sig2_celltype'.r.bw'
	### rm bedgraph
	rm $sig2_celltype'.bedgraph'
done < raw_sig_list.txt

ls *.info.txt > info_list.txt

time Rscript /storage/home/gzx103/group/software/signorm/bin/plot_sf_fsip.R info_list.txt $mark

### check
time paste *.fisher_p.txt.upperlim.txt > all.signal.txt
ls *.fisher_p.txt.upperlim.txt | awk -F '.' '{print $1"."$2}' > signal_list.txt
time Rscript /storage/home/gzx103/scratch/vision/5end/pknorm_16lim/h3k9me3_r_log/plot_ct_indexset_violin.R all.signal.txt signal_list.txt all.signal.od.violin.png

time paste *.pknorm.txt > all.signal.pkn.txt
time Rscript /storage/home/gzx103/scratch/vision/5end/pknorm_16lim/h3k9me3_r_log/plot_ct_indexset_violin.R all.signal.pkn.txt signal_list.txt all.signal.pkn.violin.png

