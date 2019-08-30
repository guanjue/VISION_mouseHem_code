##################
### cluster union merged peaks
for i in {1..20}
do
	echo $i
	rm 'cluster_'$i'_merge.bed'
	for file in $(cat 'tf_JI_plots/cluster_'$i'_list.txt')
	do
		echo $file
		cat ../bed_files_mm10/$file'.sort.bed' >> 'cluster_'$i'_merge.bed'
	done
	sort -k1,1 -k2,2n 'cluster_'$i'_merge.bed' > 'cluster_'$i'_merge.sort.bed'
	bedtools merge -i 'cluster_'$i'_merge.sort.bed' > 'cluster_'$i'_merged.sort.bed'
	rm 'cluster_'$i'_merge.sort.bed'
	rm 'cluster_'$i'_merge.bed'
done


##################
### cluster intersect peaks 2
for i in {1..20}
do
	echo $i
	file1=$(head -1 'tf_JI_plots/cluster_'$i'_list.txt')
	bedtools intersect -a 'cluster_'$i'_merged.sort.bed' -b ../bed_files_mm10/$file1'.sort.bed' -wa -c > 'cluster_'$i'_intersect.bed'
	for file in $(tail -n+2 'tf_JI_plots/cluster_'$i'_list.txt')
	do
		echo $file
		bedtools intersect -a 'cluster_'$i'_intersect.bed' -b ../bed_files_mm10/$file'.sort.bed' > 'cluster_'$i'_intersect.bed.tmp' 
		cut -f4 'cluster_'$i'_intersect.bed.tmp' > 'cluster_'$i'_intersect.bed.tmp.count'
		paste 'cluster_'$i'_intersect.bed' 'cluster_'$i'_intersect.bed.tmp.count' > 'cluster_'$i'_intersect.bed.tmp' && mv 'cluster_'$i'_intersect.bed.tmp' 'cluster_'$i'_intersect.bed'
	done
	sort -k1,1 -k2,2n 'cluster_'$i'_intersect.bed' | awk -F '\t' -v OFS='\t' '{if ($4+$5+$6+$7+$8+$9+$10+$11+$12+$13+$14+$15+$16+$17+$18+$19+$20+$21+$22+$23 >= 2) print $1, $2, $3}' > 'cluster_'$i'_intersect.sort.bed'
	rm 'cluster_'$i'_intersect.bed'
	rm 'cluster_'$i'_intersect.bed.tmp.count'
done


##################
### cluster intersect peaks
for i in {1..20}
do
	echo $i
	file1=$(head -1 'tf_JI_plots/cluster_'$i'_list.txt')
	bedtools intersect -a 'cluster_'$i'_merged.sort.bed' -b ../bed_files_mm10/$file1'.sort.bed' -wa > 'cluster_'$i'_intersect.bed'
	for file in $(tail -n+2 'tf_JI_plots/cluster_'$i'_list.txt')
	do
		echo $file
		bedtools intersect -a 'cluster_'$i'_intersect.bed' -b ../bed_files_mm10/$file'.sort.bed' > 'cluster_'$i'_intersect.bed.tmp' && mv 'cluster_'$i'_intersect.bed.tmp' 'cluster_'$i'_intersect.bed'
	done
	sort -k1,1 -k2,2n 'cluster_'$i'_intersect.bed' > 'cluster_'$i'_intersect.sort.bed'
	rm 'cluster_'$i'_intersect.bed'
done

