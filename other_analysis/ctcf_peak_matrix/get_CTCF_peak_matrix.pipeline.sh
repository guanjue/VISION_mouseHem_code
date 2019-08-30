cd ~/group/projects/vision/CTCF_pk_moredata

### get bedgraph
qsub get_bedgraph.sh

### S3norm
qsub run_pipeline.sh

###### get bigwig
### S3norm NBP
ls S3norm_NBP_bedgraph/*.bedgraph.NBP.s3norm.bedgraph > S3norm_NBP.list.txt
qsub get_bw.S3norm_NBP.sh
### NBP
ls NBP_bedgraph/*.s3norm.NB.neglog10p.bedgraph > NBP.list.txt
qsub get_bw.NBP.sh
### S3norm rc
ls S3norm_rc_bedgraph/*.bedgraph.s3norm.bedgraph > S3norm_rc.list.txt
qsub get_bw.S3norm_rc.sh

### get peaks
cat chr11_loci.txt | awk -F ' ' -v OFS='\t' '{print "chr11",$1,$2,"chr11_"$1"_"$2}' | sort -k1,1 -k2,2n > chr11_loci.bed
### S3norm NBP
sort -k4,4 chr11_loci.bed > chr11_loci.idsort.mat.S3norm_NBP.txt
for file in $(cat S3norm_rc.list.txt)
do
	echo $file
	time ~/group/software/ucsc/bigWigAverageOverBed $file'.bw' chr11_loci.bed $file'.tab'
	sort -k1,1 $file'.tab' | cut -f6 > $file'.tab.txt'
	paste chr11_loci.idsort.mat.S3norm_NBP.txt $file'.tab.txt' > chr11_loci.idsort.mat.S3norm_NBP.txt.tmp && mv chr11_loci.idsort.mat.S3norm_NBP.txt.tmp chr11_loci.idsort.mat.S3norm_NBP.txt
	rm $file'.tab.txt'
	rm $file'.tab'
done
### NBP
sort -k4,4 chr11_loci.bed > chr11_loci.idsort.mat.NBP.txt
for file in $(cat NBP.list.txt)
do
	echo $file
	time ~/group/software/ucsc/bigWigAverageOverBed $file'.bw' chr11_loci.bed $file'.tab'
	sort -k1,1 $file'.tab' | cut -f6 > $file'.tab.txt'
	paste chr11_loci.idsort.mat.NBP.txt $file'.tab.txt' > chr11_loci.idsort.mat.NBP.txt.tmp && mv chr11_loci.idsort.mat.NBP.txt.tmp chr11_loci.idsort.mat.NBP.txt
	rm $file'.tab.txt'
	rm $file'.tab'
done
### S3norm rc
sort -k4,4 chr11_loci.bed > chr11_loci.idsort.mat.S3norm_rc.txt
for file in $(cat S3norm_rc.list.txt)
do
	echo $file
	time ~/group/software/ucsc/bigWigAverageOverBed $file'.bw' chr11_loci.bed $file'.tab'
	sort -k1,1 $file'.tab' | cut -f6 > $file'.tab.txt'
	paste chr11_loci.idsort.mat.S3norm_rc.txt $file'.tab.txt' > chr11_loci.idsort.mat.S3norm_rc.txt.tmp && mv chr11_loci.idsort.mat.S3norm_rc.txt.tmp chr11_loci.idsort.mat.S3norm_rc.txt
	rm $file'.tab.txt'
	rm $file'.tab'
done


