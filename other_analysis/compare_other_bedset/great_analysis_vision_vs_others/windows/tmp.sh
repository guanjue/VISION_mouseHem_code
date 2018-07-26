echo all 4 sets
mkdir all_4
cp *.01.bed all_4/
cd all_4
rm *X*.01.bed
ls *.01.bed > all_4_list.txt

for ep300_t in $(cat /Users/universe/Documents/2018_BG/great_analysis_vision_vs_others/windows/EP300_list.txt)
do
	ct=$(echo "$ep300_t" | awk -F '.' '{print $2}')
	echo $ct
	for pattern in $(cat all_4_list.txt)
	do
		pattern_4=$(echo "$pattern" | awk -F '.' '{print $1}')
		bedtools window -a '/Users/universe/Documents/2018_BG/great_analysis_vision_vs_others/windows/'$ep300_t -b $pattern_4'.01.bed' -w 150 -u > $pattern_4'.'$ct'.bed'
		Rscript /Users/universe/Documents/2018_BG/great_analysis_vision_vs_others/windows/expect_vs_obs.R $pattern_4'.'$ct'.bed' $pattern_4'.01.bed' '/Users/universe/Documents/2018_BG/great_analysis_vision_vs_others/windows/'$ep300_t 2725521370 $pattern_4'.'$ct'.all4'
	done
	mkdir 'EP300_all4_'$ct
	mv *'_'*'_'*'_'*'.'$ct'.bed' 'EP300_all4_'$ct
done

mkdir all4_enrich
mv *.all4.enrichment.txt all4_enrich

echo get heatmap
cd all4_enrich
ls *all4.enrichment.txt | awk -F '.' '{print $1}' | sort -u > intersect_pattern.txt
for pattern in $(cat intersect_pattern.txt)
do
	echo $pattern
	tail -n+2 $pattern'.CH12.all4.enrichment.txt' >> CH12.enrichment.txt
	tail -n+2 $pattern'.Fetal_Liver.all4.enrichment.txt' >> Fetal_Liver.enrichment.txt
	tail -n+2 $pattern'.MEL.all4.enrichment.txt' >> MEL.enrichment.txt
done
paste intersect_pattern.txt CH12.enrichment.txt Fetal_Liver.enrichment.txt MEL.enrichment.txt > all4.enrich.matrix.txt
rm CH12.enrichment.txt Fetal_Liver.enrichment.txt MEL.enrichment.txt

echo plot heatmap
Rscript /Users/universe/Documents/2018_BG/great_analysis_vision_vs_others/windows/get_enrichment_heatmap.R all4.enrich.matrix.txt all4.enrich.matrix.png
