###### download data
mkdir bw_files
while read LINE
do
	file1=$(echo "$LINE" | awk -F '\t' -v OFS='\t' '{print $1}')
	file2=$(echo "$LINE" | awk -F '\t' -v OFS='\t' '{print $2}')
	mk=$(echo "$LINE" | awk -F '\t' -v OFS='\t' '{print $3}')
	echo $file1
	echo $file2
	echo $mk
	wget http://usevision.org/data/IDEASmouseHem/raw-signal/$file1
	wget http://usevision.org/data/IDEASmouseHem/raw-signal/$file2
	mv $file1 bw_files/
	mv $file2 bw_files/
done < ER4_bw_list.txt

###### get bigwigaverageoverbed
cat /storage/home/gzx103/group/projects/3d/G1E-ER4-uninduced.merged/10kb/TAD_regions.txt | awk -F ' ' -v OFS='\t' '{print $1, $2, $3, "chr"$4}' > TAD_regions.bed
cat TAD_regions.bed | awk -F '\t' -v OFS=' ' '{print $1, $2, $3, $4}' > TAD_regions.bed.txt
bedtools intersect -a /storage/home/gzx103/group/projects/vision/mm10.10KB.noblacklist.bed -b TAD_regions.bed -v | sort -k1,1 -k2,2n > mm10.10KB.noblacklist.notad.bg.bed
#cat /storage/home/gzx103/group/projects/vision/mm10.10KB.noblacklist.bed | sort -k1,1 -k2,2n > mm10.10KB.noblacklist.bg.bed
Rscript get_random_lines.R mm10.10KB.noblacklist.bg.bed bg.bed 2018 10000
bed_file=TAD_regions.bed
bed_file_bg=bg.bed
script_dir='/storage/home/g/gzx103/group/software/vision_pipeline_mouse_mm10/src/'

cp /storage/home/gzx103/group/projects/vision/merged_input/200_noblack.5bins_bgsig_mean.round.bw bw_files/
### get rc
~/group/software/ucsc/bigWigAverageOverBed bw_files/200_noblack.5bins_bgsig_mean.round.bw $bed_file bw_files/input_sig.tab
### sort bigWigAverageOverBed output
cat bw_files/input_sig.tab | awk -F '_' -v OFS='\t' '{print $1,$2,$3}' | sort -k1,1 -k2,2n | cut -f8 > bw_files/input_sig.tab.txt

### get NBP
while read LINE
do
	sig1=$(echo "$LINE" | awk -F '\t' -v OFS='\t' '{print $1}')
	sig2=$(echo "$LINE" | awk -F '\t' -v OFS='\t' '{print $2}')
	mk=$(echo "$LINE" | awk -F '\t' -v OFS='\t' '{print $3}')
	echo $sig1
	echo $sig2
	paste bw_files/$mk'.rep1.tab.txt' bw_files/$mk'.rep2.tab.txt' | awk -F '\t' -v OFS='\t' '{print ($1+$2)/2}' > bw_files/$mk'.mean_rc.txt'
done < ER4_bw_list.txt


### get NBP
mkdir nbp
while read LINE
do
	sig1=$(echo "$LINE" | awk -F '\t' -v OFS='\t' '{print $1}')
	sig2=$(echo "$LINE" | awk -F '\t' -v OFS='\t' '{print $2}')
	mk=$(echo "$LINE" | awk -F '\t' -v OFS='\t' '{print $3}')
	echo $sig1
	echo $sig2
	### get rc
	#~/group/software/ucsc/bigWigAverageOverBed 'bw_files/'$sig1 $bed_file 'bw_files/'$mk'.rep1.tab'
	#~/group/software/ucsc/bigWigAverageOverBed 'bw_files/'$sig2 $bed_file 'bw_files/'$mk'.rep2.tab'
	### sort bigWigAverageOverBed output
	#cat 'bw_files/'$mk'.rep1.tab' | awk -F '_' -v OFS='\t' '{print $1,$2,$3}' | sort -k1,1 -k2,2n | cut -f8 > 'bw_files/'$mk'.rep1.tab.txt'
	#cat 'bw_files/'$mk'.rep2.tab' | awk -F '_' -v OFS='\t' '{print $1,$2,$3}' | sort -k1,1 -k2,2n | cut -f8 > 'bw_files/'$mk'.rep2.tab.txt'
	### get rc bg
	#~/group/software/ucsc/bigWigAverageOverBed 'bw_files/'$sig1 $bed_file_bg 'bw_files/'$mk'.rep1.bg.tab'
	#~/group/software/ucsc/bigWigAverageOverBed 'bw_files/'$sig2 $bed_file_bg 'bw_files/'$mk'.rep2.bg.tab'
	### sort bigWigAverageOverBed output bg
	#cat 'bw_files/'$mk'.rep1.bg.tab' | awk -F '_' -v OFS='\t' '{print $1,$2,$3}' | sort -k1,1 -k2,2n | cut -f8 > 'bw_files/'$mk'.rep1.tab.bg.txt'
	#cat 'bw_files/'$mk'.rep2.bg.tab' | awk -F '_' -v OFS='\t' '{print $1,$2,$3}' | sort -k1,1 -k2,2n | cut -f8 > 'bw_files/'$mk'.rep2.tab.bg.txt'
	### get nbp
	Rscript 'negative_binomial_p_2r_bgadj_bayes.R' $mk'.rep1.tab.txt' 'bw_files/' input_sig.tab.txt 'bw_files/' $mk'.rep1.nbp.txt' 'bw_files/'$mk'.rep1.tab.bg.txt'
	Rscript 'negative_binomial_p_2r_bgadj_bayes.R' $mk'.rep2.tab.txt' 'bw_files/' input_sig.tab.txt 'bw_files/' $mk'.rep2.nbp.txt' 'bw_files/'$mk'.rep2.tab.bg.txt'
	mv $mk'.rep1.nbp.txt.nbp_1r_bgadj.txt' nbp/
	mv $mk'.rep2.nbp.txt.nbp_1r_bgadj.txt' nbp/
done < ER4_bw_list.txt


### merge NBP
mkdir merge_nbp
while read LINE
do
	sig1=$(echo "$LINE" | awk -F '\t' -v OFS='\t' '{print $1}')
	sig2=$(echo "$LINE" | awk -F '\t' -v OFS='\t' '{print $2}')
	mk=$(echo "$LINE" | awk -F '\t' -v OFS='\t' '{print $3}')
	echo $sig1
	echo $sig2
	Rscript $script_dir'fisher_pval.R' $mk '.nbp_1r_bgadj.txt' 'nbp/' 100
	mv $mk'.fisher_p.frip_snr.txt' merge_nbp/
	mv $mk'.fisher_p.txt' merge_nbp/
done < ER4_bw_list.txt
head merge_nbp/*.fisher_p.frip_snr.txt

### S3norm
while read LINE
do
	sig1=$(echo "$LINE" | awk -F '\t' -v OFS='\t' '{print $1}')
	sig2=$(echo "$LINE" | awk -F '\t' -v OFS='\t' '{print $2}')
	mk=$(echo "$LINE" | awk -F '\t' -v OFS='\t' '{print $3}')
	upperlim=100
	lowerlim=0
	echo $sig1 
	echo $sig2
	echo $mk
	### set upper limit
	cat 'merge_nbp/'$mk'.fisher_p.txt' | awk -F '\t' -v OFS='\t' -v ul=$upperlim '{if ($1>=ul) print ul; else print $1}' > 'merge_nbp/'$mk'.fisher_p.lim.txt'
	cat merge_nbp/h3k36me3.fisher_p.txt | awk -F '\t' -v OFS='\t' -v ul=$upperlim '{if ($1>=ul) print ul; else print $1}' > merge_nbp/ref.fisher_p.lim.txt
	### peak norm
	time python 'peaknorm_rotate_log_ref_mean.py' -n 7619 -a merge_nbp/ref.fisher_p.lim.txt -b 'merge_nbp/'$mk'.fisher_p.lim.txt' -u $upperlim -l 0
	### set lim
	cat 'merge_nbp/'$mk'.fisher_p.pknorm.ref.txt' | awk -F '\t' -v OFS='\t' -v ul=16 -v ll=3 '{if ($1>=ul) print ul; else if ($1<=ll) print ll; else print $1}' > 'merge_nbp/'$mk'.fisher_p.pknorm.ref.lim.txt'
	### rm tmp files
	#rm 'merge_nbp/'$mk'.fisher_p.lim.txt'
	#rm merge_nbp/ref.fisher_p.lim.txt
done < ER4_bw_list.txt


time bash run_IDEAS_tad.sh

