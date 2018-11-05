cd /Users/universe/Documents/2018_BG/vision_mouse_analysis/ccRE_intersect_number

mkdir bed_files

### get master peak list
rm all_pks_concate.bed
while read LINE
do
	file=$(echo "$LINE" | awk -F '\t' '{print $1}')
	echo $file
	cut -f1,2,3 $file > $file'.bed'
	cat $file'.bed' | awk -F '\t' -v OFS='\t' '{print $1, int(($2+$3)/2), int(($2+$3)/2)+1}' > $file'.mid.bed'
	cat $file'.bed' >> all_pks_concate.bed
	mv $file'.bed' bed_files
	mv $file'.mid.bed' bed_files
done < bed_list_all.txt


### merge pk bed files
sort -k1,1 -k2,2n -u all_pks_concate.bed | tail -n+2 > all_pks_concate.sort.bed
cat all_pks_concate.sort.bed | awk -F '\t' -v OFS='\t' '{print $1, int(($2+$3)/2), int(($2+$3)/2)+1}' | sort -k1,1 -k2,2n > all_pks_concate.sort.mid.bed
bedtools merge -i all_pks_concate.sort.mid.bed -d 200 > all_pks_merge.sort.bed
#cp all_pks_concate.sort.bed all_pks_merge.sort.bed

### count mat
cat all_pks_merge.sort.bed | awk -F '\t' -v OFS='\t' '{print $1"_"$2"_"$3}' > allpk_count_mat.txt

while read LINE
do
	file=$(echo "$LINE" | awk -F '\t' '{print $1}')
	echo $file
	bedtools window -a all_pks_merge.sort.bed -b 'bed_files/'$file'.mid.bed' -c -w 200 > 'bed_files/'$file'.count.bed'
	cut -f4 'bed_files/'$file'.count.bed' > 'bed_files/'$file'.count.txt'
	paste allpk_count_mat.txt 'bed_files/'$file'.count.txt' > allpk_count_mat.txt.tmp && mv allpk_count_mat.txt.tmp allpk_count_mat.txt
done < bed_list_all.txt


time Rscript /Users/universe/Documents/2018_BG/vision_pipeline_mouse_mm10/other_analysis/ccRE_intersect_number/upset.R

