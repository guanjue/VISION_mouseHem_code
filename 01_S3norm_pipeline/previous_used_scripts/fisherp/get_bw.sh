#PBS -l nodes=2:ppn=8
#PBS -l walltime=20:00:00
#PBS -j oe
#PBS -A open

module load gcc
module load python/2.7.8

cd /storage/home/gzx103/scratch/vision/5end/fisher_p_100lim

#time sleep 20000
### cp p-vlaue files
#for mark in $(cat /storage/home/gzx103/group/projects/vision/mark_list.txt)
#do
#        cp '/storage/home/gzx103/scratch/vision/5end/'$mark'/'*'.bamtobed5endintersect.signal.nbp_bgadj.txt.nbp_2r_bgadj.txt' ./
#done

### get sample list
ls *.bamtobed5endintersect.signal.nbp_bgadj.txt.nbp_2r_bgadj.txt > sample_list.txt
### extrac cell type marker list
cat sample_list.txt | awk -F '.' -v OFS='\t' '{print $1"."$2}' | sort -u > cell_marker_list.txt

### get Fisher's method combined pval
for cm in $(cat cell_marker_list.txt)
do
        echo $cm
	#Rscript /storage/home/gzx103/group/software/signorm/bin/fisher_pval.R $cm '.bamtobed5endintersect.signal.nbp_bgadj.txt.nbp_2r_bgadj.txt' '/storage/home/gzx103/scratch/vision/5end/fisher_p/' 100
done

#qsub signorm_atac.sh
#qsub signorm_ctcf.sh
#qsub signorm_h3k4me1.sh
#qsub signorm_h3k4me3.sh
#qsub signorm_h3k9me3.sh
#qsub signorm_h3k27me3.sh
#qsub signorm_h3k27ac.sh
#qsub signorm_h3k36me3.sh

### extract Fisher's pval list
ls *.fisher_p.txt > fisher_p_list.txt
### convert sig to bigwig
for filename in $(cat fisher_p_list.txt)
do
        echo $filename
        ### get bedgraph
        paste ~/group/projects/vision/merged_input/200_noblack.11_22_2017.bed $filename | awk -F '\t' -v OFS='\t' '{print $1,$2,$3,$4}' | sort -k1,1 -k2,2n> $filename'.fisher_p.bedgraph'
        ### bedGraphToBigWig
        /storage/home/gzx103/group/projects/vision/input_norm/bedGraphToBigWig $filename'.fisher_p.bedgraph' /storage/home/gzx103/group/projects/vision/input_norm/mm10.chrom.sizes $filename'.fisher_p.bw'
        ### rm bedgraph
        rm $filename'.fisher_p.bedgraph'
done

mkdir fisher_p_bw_100lim
mv *.fisher_p.bw fisher_p_bw_100lim/

#mkdir /storage/home/gzx103/scratch/vision/fisher_p_bw
#mv *.fisher_p.bw /storage/home/gzx103/scratch/vision/fisher_p_bw/

#qsub normalization_across_marks.sh

#cd /storage/home/gzx103/group/projects/vision/fisher_p_signorm/
#bash PBS_sub_fisherp.sh
