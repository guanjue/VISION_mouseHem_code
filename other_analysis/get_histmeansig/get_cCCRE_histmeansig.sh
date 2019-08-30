cd ~/group/projects/vision/get_cCER_histmeansig

cut -f1,2,3 ~/group/projects/vision/index_caller/atac_20cell.index.matrix.no0.txt | awk -F '\t' -v OFS='\t' '{print $1,$2,$3,$1"_"$2"_"$3}' | sort -k1,1 > atac_pk_no0.bed


declare -a mk_list=("H3K27ac" "H3K4me1" "H3K4me3")

for mk in "${mk_list[@]}"
do
echo $mk
cp atac_pk_no0.bed 'atac_pk_no0.bed.tab.'$mk'.txt'
for file in $(cat $mk'.list.txt')
do
	echo $file
	time ~/group/software/ucsc/bigWigAverageOverBed ~/group/projects/vision/bw/$file atac_pk_no0.bed atac_pk_no0.bed.tab.tmp
	sort -k1,1 atac_pk_no0.bed.tab.tmp | cut -f5 > atac_pk_no0.bed.tab.tmp.txt
	paste 'atac_pk_no0.bed.tab.'$mk'.txt' atac_pk_no0.bed.tab.tmp.txt > 'atac_pk_no0.bed.tab.'$mk'.txt.tmp'
	mv 'atac_pk_no0.bed.tab.'$mk'.txt.tmp' 'atac_pk_no0.bed.tab.'$mk'.txt'
	sort -k1,1 -k2,2n 'atac_pk_no0.bed.tab.'$mk'.txt' > 'atac_pk_no0.bed.tab.'$mk'.sort.txt'
done
done


