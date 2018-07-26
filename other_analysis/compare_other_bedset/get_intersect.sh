sort -k1,1 -k2,2n AmitFriedman_BloodCellEnhancerCatalog_2014_mm10.bed.txt > A.sort.bed
sort -k1,1 -k2,2n Known_erythroid_REs_2017_mm10.txt > K.sort.bed
sort -k1,1 -k2,2n atac_20cell.fun.no0.bed.txt > V.sort.bed
sort -k1,1 -k2,2n FLD_PEAKS_merged_mm10.txt > F.sort.bed
sort -k1,1 -k2,2n SCREEN_mouse_ccREs_mm10.txt > S.sort.bed
sort -k1,1 -k2,2n EP300_CH12_rep1intsxrep2.txt > EP300.CH12.sort.bed
sort -k1,1 -k2,2n EP300_MEL_rep1intsxrep2.txt > EP300.MEL.sort.bed
sort -k1,1 -k2,2n EP300_mFL14_5_rep1intsxrep2.txt > EP300.mFL14_5.sort.bed

### merge EP300
cat EP300.CH12.sort.bed EP300.MEL.sort.bed EP300.mFL14_5.sort.bed | sort -k1,1 -k2,2n > EP300_cat.sort.bed
bedtools merge -i EP300_cat.sort.bed > EP300_merged.sort.bed

### VAS_EP300
bedtools intersect -a V.sort.bed -b A.sort.bed -wa -u > VA.bed
bedtools intersect -a VA.bed -b S.sort.bed -wa -u > VAS.bed
bedtools intersect -a VAS.bed -b EP300_merged.sort.bed -wa -u > VAS_EP300.bed

### V_EP300
bedtools intersect -a V.sort.bed -b EP300_merged.sort.bed -wa -u > V_EP300.bed
bedtools intersect -a V.sort.bed -b EP300_CH12.sort.bed -wa -u > V_EP300_CH12.bed
bedtools intersect -a V.sort.bed -b EP300_MEL.sort.bed -wa -u > V_EP300_MEL.bed
bedtools intersect -a V.sort.bed -b EP300_mFL14_5.sort.bed -wa -u > V_EP300_mFL14_5.bed
wc -l V_EP300*.bed

### A_EP300, F_EP300, S_EP300
bedtools intersect -a A.sort.bed -b EP300_merged.sort.bed -wa -u > A_EP300.bed
bedtools intersect -a F.sort.bed -b EP300_merged.sort.bed -wa -u > F_EP300.bed
bedtools intersect -a S.sort.bed -b EP300_merged.sort.bed -wa -u > S_EP300.bed
wc -l *_EP300.bed

### A_K, F_K, S_K, V_K
bedtools intersect -a A.sort.bed -b K.sort.bed -wa -u > A_K.bed
bedtools intersect -a F.sort.bed -b K.sort.bed -wa -u > F_K.bed
bedtools intersect -a S.sort.bed -b K.sort.bed -wa -u > S_K.bed
bedtools intersect -a V.sort.bed -b K.sort.bed -wa -u > V_K.bed
wc -l *_K.bed



echo V_F_S

echo level uniq
echo V uniq
bedtools intersect -a V.sort.bed -b F.sort.bed -v > V_X_0_X.bed
bedtools intersect -a V_X_0_X.bed -b S.sort.bed -v > V_X_0_0.bed
echo F uniq
bedtools intersect -a F.sort.bed -b V.sort.bed -v > 0_X_F_X.bed
bedtools intersect -a 0_X_F_X.bed -b S.sort.bed -v > 0_X_F_0.bed
echo S uniq
bedtools intersect -a S.sort.bed -b V.sort.bed -v > 0_X_X_S.bed
bedtools intersect -a 0_X_X_S.bed -b F.sort.bed -v > 0_X_0_S.bed


echo level 2 uniq
echo V_X_F_0 
bedtools intersect -a V.sort.bed -b F.sort.bed -u > V_X_F_X.bed
bedtools intersect -a V_X_F_X.bed -b S.sort.bed -v > V_X_F_0.bed
echo V_X_0_S 
bedtools intersect -a V_X_0_X.bed -b S.sort.bed -u > V_X_0_S.bed
echo 0_X_F_S
bedtools intersect -a F.sort.bed -b S.sort.bed -u > X_X_F_S.bed
bedtools intersect -a X_X_F_S.bed -b V.sort.bed -v > 0_X_F_S.bed

echo level 3 uniq
echo V_X_F_X
bedtools intersect -a V_X_F_X.bed -b S.sort.bed -u > V_X_F_S.bed



cut -f1,2,3 V_X_0_0.bed > V_X_0_0.01.bed
cut -f1,2,3 0_X_F_0.bed > 0_X_F_0.01.bed
cut -f1,2,3 0_X_0_S.bed > 0_X_0_S.01.bed
cut -f1,2,3 V_X_F_0.bed > V_X_F_0.01.bed
cut -f1,2,3 V_X_0_S.bed > V_X_0_S.01.bed
cut -f1,2,3 0_X_F_S.bed > 0_X_F_S.01.bed
cut -f1,2,3 V_X_F_S.bed > V_X_F_S.01.bed

echo get intersect with known ccREs
bedtools intersect -a Known_erythroid_REs_2017_mm10.txt -b V_X_0_0.01.bed -u > V_X_0_0.01.K.bed
bedtools intersect -a Known_erythroid_REs_2017_mm10.txt -b 0_X_F_0.01.bed -u > 0_X_F_0.01.K.bed
bedtools intersect -a Known_erythroid_REs_2017_mm10.txt -b 0_X_0_S.01.bed -u > 0_X_0_S.01.K.bed
bedtools intersect -a Known_erythroid_REs_2017_mm10.txt -b V_X_F_0.01.bed -u > V_X_F_0.01.K.bed
bedtools intersect -a Known_erythroid_REs_2017_mm10.txt -b V_X_0_S.01.bed -u > V_X_0_S.01.K.bed
bedtools intersect -a Known_erythroid_REs_2017_mm10.txt -b 0_X_F_S.01.bed -u > 0_X_F_S.01.K.bed
bedtools intersect -a Known_erythroid_REs_2017_mm10.txt -b V_X_F_S.01.bed -u > V_X_F_S.01.K.bed

wc -l *_X_*_*.01.K.bed

ls *_X_*_*.01.bed > VFS_list.txt
#ls *_X_*_*.01.K.bed | awk -F '.' '{print $1}' | awk -F '_' -v OFS='_' '{print $1,$3,$4, "K"}' >> VFS_list.txt

echo get ven diagram
python plot_ven.py -l VFS_list.txt

for file in $(cat VFS_list.txt)
do
	echo $file
	pattern=$(echo "$file" | awk -F '.' '{print $1}')
	bedtools intersect -a Known_erythroid_REs_2017_mm10.txt -b $pattern'.01.bed' -u > $pattern'.01.K.bed'
done


echo intersect with EP300
echo Download bed narrow peaks (conservative idr thresholded peaks) from ENCODE
echo unzip *.gz files

ls EP300.*.sort.bed > EP300_list.txt

for ep300_t in $(cat EP300_list.txt)
do
	ct=$(echo "$ep300_t" | awk -F '.' '{print $2}')
	echo $ct
	bedtools intersect -a $ep300_t -b V_X_0_0.01.bed -u > 'V_X_0_0.'$ct'.bed'
	Rscript expect_vs_obs.R 'V_X_0_0.'$ct'.bed' V_X_0_0.01.bed $ep300_t 2725521370 'V_X_0_0.'$ct'.VFS'
	bedtools intersect -a $ep300_t -b 0_X_F_0.01.bed -u > '0_X_F_0.'$ct'.bed'
	Rscript expect_vs_obs.R '0_X_F_0.'$ct'.bed' 0_X_F_0.01.bed $ep300_t 2725521370 '0_X_F_0.'$ct'.VFS'
	bedtools intersect -a $ep300_t -b 0_X_0_S.01.bed -u > '0_X_0_S.'$ct'.bed'
	Rscript expect_vs_obs.R '0_X_0_S.'$ct'.bed' 0_X_0_S.01.bed $ep300_t 2725521370 '0_X_0_S.'$ct'.VFS'
	bedtools intersect -a $ep300_t -b V_X_F_0.01.bed -u > 'V_X_F_0.'$ct'.bed'
	Rscript expect_vs_obs.R 'V_X_F_0.'$ct'.bed' V_X_F_0.01.bed $ep300_t 2725521370 'V_X_F_0.'$ct'.VFS'
	bedtools intersect -a $ep300_t -b V_X_0_S.01.bed -u > 'V_X_0_S.'$ct'.bed'
	Rscript expect_vs_obs.R 'V_X_0_S.'$ct'.bed' V_X_0_S.01.bed $ep300_t 2725521370 'V_X_0_S.'$ct'.VFS'
	bedtools intersect -a $ep300_t -b 0_X_F_S.01.bed -u > '0_X_F_S.'$ct'.bed'
	Rscript expect_vs_obs.R '0_X_F_S.'$ct'.bed' 0_X_F_S.01.bed $ep300_t 2725521370 '0_X_F_S.'$ct'.VFS'
	bedtools intersect -a $ep300_t -b V_X_F_S.01.bed -u > 'V_X_F_S.'$ct'.bed'
	Rscript expect_vs_obs.R 'V_X_F_S.'$ct'.bed' V_X_F_S.01.bed $ep300_t 2725521370 'V_X_F_S.'$ct'.VFS'
	mkdir 'EP300_VFS_'$ct
	mv *'_'*'_'*'_'*'.'$ct'.bed' 'EP300_VFS_'$ct
done

mkdir VFS_enrich
mv *.VFS.enrichment.txt VFS_enrich

echo get heatmap
cd VFS_enrich
ls *VFS.enrichment.txt | awk -F '.' '{print $1}' | sort -u > intersect_pattern.txt
for pattern in $(cat intersect_pattern.txt)
do
	echo $pattern
	tail -n+2 $pattern'.CH12.VFS.enrichment.txt' >> CH12.enrichment.txt
	tail -n+2 $pattern'.mFL14_5.VFS.enrichment.txt' >> mFL14_5.enrichment.txt
	tail -n+2 $pattern'.MEL.VFS.enrichment.txt' >> MEL.enrichment.txt
done


paste intersect_pattern.txt CH12.enrichment.txt mFL14_5.enrichment.txt MEL.enrichment.txt > VFS.enrich.matrix.txt
rm CH12.enrichment.txt mFL14_5.enrichment.txt MEL.enrichment.txt

echo plot heatmap
Rscript /Users/universe/Documents/2018_BG/vision_pipeline_mouse_mm10/other_analysis/compare_other_bedset/get_enrichment_heatmap.R VFS.enrich.matrix.txt VFS.enrich.matrix.png

cd ..
echo all 4 sets
mkdir all_4
cp *.01.bed all_4/
cd all_4
#rm *X*.01.bed
ls *.01.bed > all_4_list.txt

for ep300_t in $(cat /Users/universe/Documents/2018_BG/great_analysis_vision_vs_others/windows/EP300_list.txt)
do
	ct=$(echo "$ep300_t" | awk -F '.' '{print $2}')
	echo $ct
	for pattern in $(cat all_4_list.txt)
	do
		pattern_4=$(echo "$pattern" | awk -F '.' '{print $1}')
		bedtools intersect -a '/Users/universe/Documents/2018_BG/great_analysis_vision_vs_others/windows/'$ep300_t -b $pattern_4'.01.bed' -u > $pattern_4'.'$ct'.bed'
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



