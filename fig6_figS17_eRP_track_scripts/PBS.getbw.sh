#PBS -l nodes=1:ppn=8
#PBS -l walltime=10:00:00
#PBS -j oe
#PBS -A yzz2_e_g_sc_default
#PBS -l pmem=16gb
module load gcc/5.3.1
module load python/2.7.14-anaconda5.0.1
module load bedtools
cd /storage/home/gzx103/group/projects/vision/pknorm_2_16lim_ref1mo_sample1mo_0424/get_eRP_tracks/

get_eRPs_8track_file='/storage/home/gzx103/group/projects/vision/pknorm_2_16lim_ref1mo_sample1mo_0424/get_eRP_tracks/get_eRPs.8track.R'
bedGraphToBigWig_file='~/group/software/ucsc/bedGraphToBigWig'
genome_size='~/group/genome/mm10/mm10.chrom.sort.sizes'
tss_coeffi='~/group/projects/vision/rna/tss_eRP_2kb.with0.txt'
dist_coeffi='~/group/projects/vision/rna/dist_eRP_2kb.with0.txt'

for ct in $(cat ct_list.txt)
do
	### get percentage
	cd $ct'_folder'
	time Rscript $get_eRPs_8track_file $ct $tss_coeffi $dist_coeffi
	### get bigwig
	for i in {1..4}
	do
		echo $i
		sort -k1,1 -k2,2n $ct'.tss.'$i'.bedgraph' | awk -F '\t' -v OFS='\t' '{if ($4!="NA") print $0}' > $ct'.tss.'$i'.sort.bedgraph'
		bedtools merge -i $ct'.tss.'$i'.sort.bedgraph' -c 4 -o mean > $ct'.tss.'$i'.sort.merge.bedgraph'
		time $bedGraphToBigWig_file $ct'.tss.'$i'.sort.merge.bedgraph' $genome_size $ct'.tss.'$i'.sort.bw'
	done
	for i in {1..4}
	do
		echo $i
		sort -k1,1 -k2,2n $ct'.dist.'$i'.bedgraph' | awk -F '\t' -v OFS='\t' '{if ($4!="NA") print $0}' > $ct'.dist.'$i'.sort.bedgraph'
		bedtools merge -i $ct'.dist.'$i'.sort.bedgraph' -c 4 -o mean > $ct'.dist.'$i'.sort.merge.bedgraph'
		time $bedGraphToBigWig_file $ct'.dist.'$i'.sort.merge.bedgraph' $genome_size $ct'.dist.'$i'.sort.bw'
	done
	cd ..
done


