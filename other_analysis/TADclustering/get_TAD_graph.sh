### get beds
tail -n+2 run_TAD_IDEAS_bd_result/run_TAD_IDEAS_bd.state | awk -F ' ' -v OFS='\t' '{print $2, $3, $4, $5}' > TAD_bd.bed
tail -n+2 run_TAD_IDEAS_result/run_TAD_IDEAS.state | awk -F ' ' -v OFS='\t' '{print $2, $3, $4, $5}' > TAD.bed

### merge TAD & boundary
cat TAD_bd.bed TAD.bed | sort -k1,1 -k2,2n > TAD_and_bd.bed
bedtools merge -i TAD_and_bd.bed > TAD_and_bd.merge.bed

### get each module pks bed
mkdir module_pk_list
for i in {1..856}
do
	echo $i
	tail -n+$i TAD_and_bd.merge.bed | head -1 > 'TAD_and_bd.'$i'.bed'
	### intersect bds
	bedtools intersect -a TAD_bd.bed -b 'TAD_and_bd.'$i'.bed' -wa > 'TAD_bd.'$i'.bed'
	cat 'TAD_bd.'$i'.bed' | awk -F '\t' -v OFS='\t' '{print $1, $2, $3, "bd_"$4}' > 'TAD_bd.'$i'.bed.tmp' && mv 'TAD_bd.'$i'.bed.tmp' 'TAD_bd.'$i'.bed'
	### intersect tads
	bedtools intersect -a TAD.bed -b 'TAD_and_bd.'$i'.bed' -wa > 'TAD.'$i'.bed'
	cat 'TAD.'$i'.bed' | awk -F '\t' -v OFS='\t' '{print $1, $2, $3, "tad_"$4}' > 'TAD.'$i'.bed.tmp' && mv 'TAD.'$i'.bed.tmp' 'TAD.'$i'.bed'
	### get module peak with ideas state
	cat 'TAD.'$i'.bed' 'TAD_bd.'$i'.bed' | sort -k1,1 -k2,2n > 'TAD_and_bd.'$i'.state.bed'
	### mv & rm files
	mv 'TAD_and_bd.'$i'.state.bed' module_pk_list/
	rm 'TAD.'$i'.bed' 'TAD_bd.'$i'.bed' 'TAD_and_bd.'$i'.bed'
done




