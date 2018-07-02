#PBS -l nodes=1:ppn=8
#PBS -l walltime=10:00:00
#PBS -j oe
#PBS -A open

module load gcc
module load python/2.7.8

##################################
cd /storage/home/gzx103/scratch/vision/5end/atac

for sig in $(cat atac_list.txt)
do
	echo $sig
	time Rscript ~/group/software/signorm/bin/negative_binomial_p_2r_bgadj.R $sig ~/group/projects/vision/merged_input/ones_input.ataconly.rounding.txt $sig'.nbp_bgadj.txt'
	echo $sig
	###### bw file of raw_5end
	### get bedgraph
	paste ~/group/projects/vision/merged_input/200_noblack.11_22_2017.bed $sig | awk -F '\t' -v OFS='\t' '{print $1,$2,$3,$4}' | sort -k1,1 -k2,2n> $sig'.bedgraph'
	### bedGraphToBigWig
	/storage/home/gzx103/group/projects/vision/input_norm/bedGraphToBigWig $sig'.bedgraph' /storage/home/gzx103/group/projects/vision/input_norm/mm10.chrom.sizes $sig'.bw'
	### rm bedgraph
	#rm $sig'.bedgraph'
	###### bw file of nbp_nobg_sample 100lim
	### get bedgraph
	paste ~/group/projects/vision/merged_input/200_noblack.11_22_2017.bed $sig'.nbp_bgadj.txt.nbp_2r_bgadj.txt' | awk -F '\t' -v OFS='\t' '{print $1,$2,$3,$4}' | sort -k1,1 -k2,2n> $sig'.nbp_bgadj.txt.nbp_2r_bgadj.bedgraph'
	### bedGraphToBigWig
	/storage/home/gzx103/group/projects/vision/input_norm/bedGraphToBigWig $sig'.nbp_bgadj.txt.nbp_2r_bgadj.bedgraph' /storage/home/gzx103/group/projects/vision/input_norm/mm10.chrom.sizes $sig'.nbp_2r_bgadj.bw'
	### rm bedgraph
	rm $sig'.nbp_bgadj.txt.nbp_2r_bgadj.bedgraph'
	###### bw file of nbp_nobg_sample 16lim
	paste ~/group/projects/vision/merged_input/200_noblack.11_22_2017.bed $sig'.nbp_bgadj.txt.nbp_2r_bgadj.txt' | awk -F '\t' -v OFS='\t' '{if ($4<=16) print $1,$2,$3,$4; else print $1,$2,$3,16}' | sort -k1,1 -k2,2n> $sig'.nbp_bgadj.txt.nbp_2r_bgadj.16lim.bedgraph'
	### bedGraphToBigWig
	/storage/home/gzx103/group/projects/vision/input_norm/bedGraphToBigWig $sig'.nbp_bgadj.txt.nbp_2r_bgadj.16lim.bedgraph' /storage/home/gzx103/group/projects/vision/input_norm/mm10.chrom.sizes $sig'.nbp_2r_bgadj.16lim.bw'
	### rm bedgraph
	rm $sig'.nbp_bgadj.txt.nbp_2r_bgadj.16lim.bedgraph'

done

