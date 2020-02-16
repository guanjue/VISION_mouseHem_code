#PBS -l nodes=1:ppn=8
#PBS -l walltime=10:00:00
#PBS -j oe
#PBS -A open

module load gcc
module load python/2.7.8

##################################
cd /storage/home/gzx103/scratch/vision/5end/h3k4me1

while read LINE
do
	sig1=$(echo "$LINE" | awk '{print $1}')
	sig2=$(echo "$LINE" | awk '{print $2}')
	echo $sig1 $sig2
	#time Rscript ~/group/software/signorm/bin/negative_binomial_p_2r_bgadj.R $sig1 $sig2 $sig1'.nbp_bgadj.txt'
        echo $sig1
        ### get bedgraph
        #paste ~/group/projects/vision/merged_input/200_noblack.11_22_2017.bed $sig1 | awk -F '\t' -v OFS='\t' '{print $1,$2,$3,$4}' | sort -k1,1 -k2,2n> $sig1'.bedgraph'
        ### bedGraphToBigWig
        #/storage/home/gzx103/group/projects/vision/input_norm/bedGraphToBigWig $sig1'.bedgraph' /storage/home/gzx103/group/projects/vision/input_norm/mm10.chrom.sizes $sig1'.bw'
        ### rm bedgraph
        #rm $sig1'.bedgraph'
        ###### bw file of nbp_nobg_sample 100lim
        ### get bedgraph
        paste ~/group/projects/vision/merged_input/200_noblack.11_22_2017.bed $sig1'.nbp_bgadj.txt.nbp_2r_bgadj.txt' | awk -F '\t' -v OFS='\t' '{print $1,$2,$3,$4}' | sort -k1,1 -k2,2n> $sig1'.nbp_bgadj.txt.nbp_2r_bgadj.bedgraph'
        ### bedGraphToBigWig
        /storage/home/gzx103/group/projects/vision/input_norm/bedGraphToBigWig $sig1'.nbp_bgadj.txt.nbp_2r_bgadj.bedgraph' /storage/home/gzx103/group/projects/vision/input_norm/mm10.chrom.sizes $sig1'.nbp_2r_bgadj.bw'
        ### rm bedgraph
        rm $sig1'.nbp_bgadj.txt.nbp_2r_bgadj.bedgraph'
        ###### bw file of nbp_nobg_sample 16lim
        paste ~/group/projects/vision/merged_input/200_noblack.11_22_2017.bed $sig1'.nbp_bgadj.txt.nbp_2r_bgadj.txt' | awk -F '\t' -v OFS='\t' '{if ($4<=16) print $1,$2,$3,$4; else print $1,$2,$3,16}' | sort -k1,1 -k2,2n> $sig1'.nbp_bgadj.txt.nbp_2r_bgadj.16lim.bedgraph'
        ### bedGraphToBigWig
        /storage/home/gzx103/group/projects/vision/input_norm/bedGraphToBigWig $sig1'.nbp_bgadj.txt.nbp_2r_bgadj.16lim.bedgraph' /storage/home/gzx103/group/projects/vision/input_norm/mm10.chrom.sizes $sig1'.nbp_2r_bgadj.16lim.bw'
        ### rm bedgraph
        rm $sig1'.nbp_bgadj.txt.nbp_2r_bgadj.16lim.bedgraph'
done < info_table_histone_input.txt

