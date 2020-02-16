
#PBS -l nodes=2:ppn=8
#PBS -l walltime=10:00:00
#PBS -j oe
#PBS -A yzz2_c_t_sc_default

module load gcc
module load python/2.7.8


mark=ref
cd '/storage/home/gzx103/scratch/vision/5end/pknorm_16lim/'$mark'_r_log_m'
cp /storage/home/gzx103/scratch/vision/5end/pknorm_16lim/atac_r_log_mm/MK_imm_ad.atacrep.fisher_p.txt ./
cp /storage/home/gzx103/scratch/vision/5end/pknorm_16lim/ctcf_r_log_mm/ERY_fl.ctcfrep.fisher_p.txt ./
cp /storage/home/gzx103/scratch/vision/5end/pknorm_16lim/h3k4me1_r_log_mm/ERY_fl.h3k4me1rep.fisher_p.txt ./
cp /storage/home/gzx103/scratch/vision/5end/pknorm_16lim/h3k4me3_r_log_mm/ERY_fl.h3k4me3rep.fisher_p.txt ./
cp /storage/home/gzx103/scratch/vision/5end/pknorm_16lim/h3k9me3_r_log_mm/ERY_fl.h3k9me3rep.fisher_p.txt ./
cp /storage/home/gzx103/scratch/vision/5end/pknorm_16lim/h3k27me3_r_log_mm/ERY_fl.h3k27me3rep.fisher_p.txt ./
cp /storage/home/gzx103/scratch/vision/5end/pknorm_16lim/h3k36me3_r_log_mm/G1E.h3k36me3rep.fisher_p.txt ./
cp /storage/home/gzx103/scratch/vision/5end/pknorm_16lim/h3k27ac_r_log_mm/MK_imm_ad.h3k27acrep.fisher_p.txt ./

#mkdir mean_ref
#cd mean_ref/
#cp /storage/home/gzx103/scratch/vision/5end/pknorm_16lim/atac_r_log/*.atacrep.fisher_p.txt ./
#cp /storage/home/gzx103/scratch/vision/5end/pknorm_16lim/ctcf_r_log/*.ctcfrep.fisher_p.txt ./
#cp /storage/home/gzx103/scratch/vision/5end/pknorm_16lim/h3k4me1_r_log/*.h3k4me1rep.fisher_p.txt ./
#cp /storage/home/gzx103/scratch/vision/5end/pknorm_16lim/h3k4me3_r_log/*.h3k4me3rep.fisher_p.txt ./
#cp /storage/home/gzx103/scratch/vision/5end/pknorm_16lim/h3k9me3_r_log/*.h3k9me3rep.fisher_p.txt ./
#cp /storage/home/gzx103/scratch/vision/5end/pknorm_16lim/h3k27me3_r_log/*.h3k27me3rep.fisher_p.txt ./
#cp /storage/home/gzx103/scratch/vision/5end/pknorm_16lim/h3k36me3_r_log/*.h3k36me3rep.fisher_p.txt ./
#cp /storage/home/gzx103/scratch/vision/5end/pknorm_16lim/h3k27ac_r_log/*.h3k27acrep.fisher_p.txt ./

#rm LSK_BM.h3k27me3rep.fisher_p.txt
#rm LSK_BM.h3k36me3rep.fisher_p.txt
#rm T_CD4_SPL.h3k27me3rep.fisher_p.txt
#rm T_CD4_SPL.h3k27acrep.fisher_p.txt
#rm T_CD8_SPL.h3k27acrep.fisher_p.txt

#for ref_mark in $(cat /storage/home/gzx103/scratch/vision/5end/pknorm_16lim/pcor_100lim_mean/mark_list.txt)
#do
#	echo $ref_mark
#	paste *$ref_mark'rep.fisher_p.txt' | awk -F '\t' -v OFS='\t' '{sum=0; for(i=1; i<=NF; i++){sum+=$i}; sum/=NF; print sum}' >  'mean.'$ref_mark'.txt'
#	mv 'mean.'$ref_mark'.txt' '/storage/home/gzx103/scratch/vision/5end/pknorm_16lim/'$mark'_r_log_m/'
#done
#cd ..
#rm *.upperlim.txt
#paste mean.*.txt | awk -F '\t' -v OFS='\t' '{sum=0; for(i=1; i<=NF; i++){sum+=$i}; sum/=NF; print sum}' > allref.mean.txt

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
	time python /storage/home/gzx103/group/software/signorm/bin/peaknorm_rotate_log_ref_mean.py -w /storage/home/gzx103/scratch/vision/5end/pknorm_16lim/200_noblack.11_22_2017.bed -p $mark'TablePeaksFiltered.txt' -n 500000 -a $sig1_col -b $sig1'.upperlim.txt' -c $sig2_col -d $sig2'.upperlim.txt' -u $upperlim -l $lowerlim
	### get bed files
	#paste /storage/home/gzx103/scratch/vision/5end/pknorm_16lim/200_noblack.11_22_2017.bed $sig2_celltype'.pknorm.txt' | awk -F '\t' -v OFS='\t' '{print $1, $2, $3, ".", $4}' > $sig2_celltype'.pknorm.bed'
	#paste ~/group/projects/vision/merged_input/200_noblack.11_22_2017.bed $sig2_celltype'.pknorm.txt' | awk -F '\t' -v OFS='\t' '{print $1,$2,$3,$4}' | sort -k1,1 -k2,2n> $sig2_celltype'.bedgraph'
	### bedGraphToBigWig
	#/storage/home/gzx103/group/projects/vision/input_norm/bedGraphToBigWig $sig2_celltype'.bedgraph' /storage/home/gzx103/group/projects/vision/input_norm/mm10.chrom.sizes $sig2_celltype'.r.bw'
	### rm bedgraph
	#rm $sig2_celltype'.bedgraph'
done < raw_sig_list.txt

ls mean*.info.txt > info_list.txt

time Rscript /storage/home/gzx103/group/software/signorm/bin/plot_sf_fsip.R info_list.txt $mark



cp ERY_fl_ctcfrep.pknorm.txt /storage/home/gzx103/scratch/vision/5end/pknorm_16lim/ctcf_r_log_mm/ERY_fl.ctcfrep.fisher_p.ref.txt
cp ERY_fl_h3k4me1rep.pknorm.txt /storage/home/gzx103/scratch/vision/5end/pknorm_16lim/h3k4me1_r_log_mm/ERY_fl.h3k4me1rep.fisher_p.ref.txt
cp ERY_fl_h3k4me3rep.pknorm.txt /storage/home/gzx103/scratch/vision/5end/pknorm_16lim/h3k4me3_r_log_mm/ERY_fl.h3k4me3rep.fisher_p.ref.txt
cp ERY_fl_h3k9me3rep.pknorm.txt /storage/home/gzx103/scratch/vision/5end/pknorm_16lim/h3k9me3_r_log_mm/ERY_fl.h3k9me3rep.fisher_p.ref.txt
cp ERY_fl_h3k27me3rep.pknorm.txt /storage/home/gzx103/scratch/vision/5end/pknorm_16lim/h3k27me3_r_log_mm/ERY_fl.h3k27me3rep.fisher_p.ref.txt
cp G1E_h3k36me3rep.pknorm.txt /storage/home/gzx103/scratch/vision/5end/pknorm_16lim/h3k36me3_r_log_mm/G1E.h3k36me3rep.fisher_p.ref.txt
cp MK_imm_ad_atacrep.pknorm.txt ~/scratch/vision/5end/pknorm_16lim/atac_r_log_mm/MK_imm_ad.atacrep.fisher_p.ref.txt
cp MK_imm_ad_h3k27acrep.pknorm.txt ~/scratch/vision/5end/pknorm_16lim/h3k27ac_r_log_mm/MK_imm_ad.h3k27acrep.fisher_p.ref.txt

#cp mean_atac.pknorm.txt /storage/home/gzx103/scratch/vision/5end/pknorm_16lim/atac_r_log_mm/mean.atacrep.fisher_p.ref.txt
#cp mean_ctcf.pknorm.txt /storage/home/gzx103/scratch/vision/5end/pknorm_16lim/ctcf_r_log_mm/mean.ctcfrep.fisher_p.ref.txt
#cp mean_h3k4me1.pknorm.txt /storage/home/gzx103/scratch/vision/5end/pknorm_16lim/h3k4me1_r_log_mm/mean.h3k4me1rep.fisher_p.ref.txt
#cp mean_h3k4me3.pknorm.txt /storage/home/gzx103/scratch/vision/5end/pknorm_16lim/h3k4me3_r_log_mm/mean.h3k4me3rep.fisher_p.ref.txt
#cp mean_h3k9me3.pknorm.txt /storage/home/gzx103/scratch/vision/5end/pknorm_16lim/h3k9me3_r_log_mm/mean.h3k9me3rep.fisher_p.ref.txt
#cp mean_h3k27me3.pknorm.txt /storage/home/gzx103/scratch/vision/5end/pknorm_16lim/h3k27me3_r_log_mm/mean.h3k27me3rep.fisher_p.ref.txt
#cp mean_h3k36me3.pknorm.txt /storage/home/gzx103/scratch/vision/5end/pknorm_16lim/h3k36me3_r_log_mm/mean.h3k36me3rep.fisher_p.ref.txt
#cp mean_h3k27ac.pknorm.txt /storage/home/gzx103/scratch/vision/5end/pknorm_16lim/h3k27ac_r_log_mm/mean.h3k27acrep.fisher_p.ref.txt

qsub /storage/home/gzx103/scratch/vision/5end/pknorm_16lim/atac_r_log_mm/pknorm_atac_r.sh
qsub /storage/home/gzx103/scratch/vision/5end/pknorm_16lim/ctcf_r_log_mm/pknorm_ctcf_r.sh
qsub /storage/home/gzx103/scratch/vision/5end/pknorm_16lim/h3k4me1_r_log_mm/pknorm_h3k4me1_r.sh
qsub /storage/home/gzx103/scratch/vision/5end/pknorm_16lim/h3k4me3_r_log_mm/pknorm_h3k4me3_r.sh
qsub /storage/home/gzx103/scratch/vision/5end/pknorm_16lim/h3k9me3_r_log_mm/pknorm_h3k9me3_r.sh
qsub /storage/home/gzx103/scratch/vision/5end/pknorm_16lim/h3k27me3_r_log_mm/pknorm_h3k27me3_r.sh
qsub /storage/home/gzx103/scratch/vision/5end/pknorm_16lim/h3k36me3_r_log_mm/pknorm_h3k36me3_r.sh
qsub /storage/home/gzx103/scratch/vision/5end/pknorm_16lim/h3k27ac_r_log_mm/pknorm_h3k27ac_r.sh


cd '/storage/home/gzx103/scratch/vision/5end/pknorm_16lim/'$mark'_r_log'
### check
time paste mean*.fisher_p.txt.upperlim.txt > all.signal.txt
ls mean*.fisher_p.txt.upperlim.txt | awk -F '.' '{print $1"."$2}' > signal_list.txt
time Rscript /storage/home/gzx103/scratch/vision/5end/pknorm_16lim/h3k9me3_r_log/plot_ct_indexset_violin.R all.signal.txt signal_list.txt all.signal.od.violin.png

time paste *.pknorm.txt > all.signal.pkn.txt
time Rscript /storage/home/gzx103/scratch/vision/5end/pknorm_16lim/h3k9me3_r_log/plot_ct_indexset_violin.R all.signal.pkn.txt signal_list.txt all.signal.pkn.violin.png

