cp ~/group/projects/vision/vision_cres.txt ./
cat vision_cres.txt | awk -F ' ' -v OFS='\t' '{print $1,$2,$3,$1"_"$2"_"$3}' > vision_cres.bed

sort -k1,1 -k2,2n vision_cres.bed > vision_cres.sort.bed
cp vision_cres.sort.bed vision_cres.mat.txt

while read LINE
do
	echo $LINE
	file=$(echo $LINE | awk -F ' ' '{print $1}')
	echo $file
	time sort -k1,1 -k2,2n $file > $file'.sort.bed'
	bedtools map -c 5 -null 0 -o mean -a vision_cres.sort.bed -b $file'.sort.bed' > tmp.count.bed
	cut -f5 tmp.count.bed > tmp.count.txt
	paste vision_cres.mat.txt tmp.count.txt > vision_cres.mat.txt.tmp && mv vision_cres.mat.txt.tmp vision_cres.mat.txt
done < signal_list.12ct.txt


chr tss gene_category strand CFU_E_ad CFUMK CMP ERY_fl GMP MK_imm_ad LSK_BM MEP MONO_BM NEU ER4 G1E

