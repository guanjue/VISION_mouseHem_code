echo V_A_F

echo level uniq
echo V uniq
bedtools window -a atac_20cell.fun.no0.bed.txt -b AmitFriedman_BloodCellEnhancerCatalog_2014_mm10.bed.txt -w 150 -v > V_0_X_X.bed
bedtools window -a V_0_X_X.bed -b FLD_PEAKS_merged_mm10.txt -w 150 -v > V_0_0_X.bed
echo A uniq
bedtools window -a AmitFriedman_BloodCellEnhancerCatalog_2014_mm10.bed.txt -b atac_20cell.fun.no0.bed.txt -w 150 -v > 0_A_X_X.bed
bedtools window -a 0_A_X_X.bed -b FLD_PEAKS_merged_mm10.txt -w 150 -v > 0_A_0_X.bed
echo F uniq
bedtools window -a FLD_PEAKS_merged_mm10.txt -b atac_20cell.fun.no0.bed.txt -w 150 -v > 0_X_F_X.bed
bedtools window -a 0_X_F_X.bed -b AmitFriedman_BloodCellEnhancerCatalog_2014_mm10.bed.txt -w 150 -v > 0_0_F_X.bed


echo level 2 uniq
echo V_A_0_X 
bedtools window -a atac_20cell.fun.no0.bed.txt -b AmitFriedman_BloodCellEnhancerCatalog_2014_mm10.bed.txt -w 150 -u > V_A_X_X.bed
bedtools window -a V_A_X_X.bed -b FLD_PEAKS_merged_mm10.txt -w 150 -v > V_A_0_X.bed
echo V_0_F_X 
bedtools window -a V_0_X_X.bed -b FLD_PEAKS_merged_mm10.txt -w 150 -u > V_0_F_X.bed
echo 0_A_F_X
bedtools window -a AmitFriedman_BloodCellEnhancerCatalog_2014_mm10.bed.txt -b FLD_PEAKS_merged_mm10.txt -w 150 -u > X_A_F_X.bed
bedtools window -a X_A_F_X.bed -b atac_20cell.fun.no0.bed.txt -w 150 -v > 0_A_F_X.bed

echo level 3 uniq
echo V_A_F_X
bedtools window -a V_A_X_X.bed -b FLD_PEAKS_merged_mm10.txt -w 150 -u > V_A_F_X.bed

echo extract bed file
cut -f1,2,3 V_0_0_X.bed > V_0_0_X.01.bed
cut -f1,2,3 0_A_0_X.bed > 0_A_0_X.01.bed
cut -f1,2,3 0_0_F_X.bed > 0_0_F_X.01.bed
cut -f1,2,3 V_A_0_X.bed > V_A_0_X.01.bed
cut -f1,2,3 V_0_F_X.bed > V_0_F_X.01.bed
cut -f1,2,3 0_A_F_X.bed > 0_A_F_X.01.bed
cut -f1,2,3 V_A_F_X.bed > V_A_F_X.01.bed

echo intersect with known ccREs
bedtools window -a Known_erythroid_REs_2017_mm10.txt -b V_0_0_X.01.bed -w 150 -u > V_0_0_X.01.K.bed
bedtools window -a Known_erythroid_REs_2017_mm10.txt -b 0_A_0_X.01.bed -w 150 -u > 0_A_0_X.01.K.bed
bedtools window -a Known_erythroid_REs_2017_mm10.txt -b 0_0_F_X.01.bed -w 150 -u > 0_0_F_X.01.K.bed
bedtools window -a Known_erythroid_REs_2017_mm10.txt -b V_A_0_X.01.bed -w 150 -u > V_A_0_X.01.K.bed
bedtools window -a Known_erythroid_REs_2017_mm10.txt -b V_0_F_X.01.bed -w 150 -u > V_0_F_X.01.K.bed
bedtools window -a Known_erythroid_REs_2017_mm10.txt -b 0_A_F_X.01.bed -w 150 -u > 0_A_F_X.01.K.bed
bedtools window -a Known_erythroid_REs_2017_mm10.txt -b V_A_F_X.01.bed -w 150 -u > V_A_F_X.01.K.bed

wc -l *_X.01.K.bed

ls *_X.01.bed > VAF_list.txt

for file in $(cat VAF_list.txt)
do
	echo $file
	pattern=$(echo "$file" | awk -F '.' '{print $1}')
	bedtools window -a Known_erythroid_REs_2017_mm10.txt -b $pattern'.01.bed' -w 150 -u > $pattern'.01.K.bed'
done


### get script for ven
wc -l *_X.01.bed > VAF_list.count.tmp.txt
head -7 VAF_list.count.tmp.txt | awk -F ' ' '{print }'


echo intersect with EP300
echo Download bed narrow peaks (conservative idr thresholded peaks) from ENCODE
echo unzip *.gz files
gunzip -k *.gz
ls ENC*.bed > EP300_list.txt
for ep300_t in $(cat EP300_list.txt)
do
	ct=$(echo "$ep300_t" | awk -F '.' '{print $2}')
	echo $ct
	bedtools window -a $ep300_t -b V_0_0_X.01.bed -w 150 -u > 'V_0_0_X.'$ct'.bed'
	Rscript expect_vs_obs.R 'V_0_0_X.'$ct'.bed' V_0_0_X.01.bed $ep300_t 2725521370 'V_0_0_X.'$ct'.VAF'
	bedtools window -a $ep300_t -b 0_A_0_X.01.bed -w 150 -u > '0_A_0_X.'$ct'.bed'
	Rscript expect_vs_obs.R '0_A_0_X.'$ct'.bed' 0_A_0_X.01.bed $ep300_t 2725521370 '0_A_0_X.'$ct'.VAF'
	bedtools window -a $ep300_t -b 0_0_F_X.01.bed -w 150 -u > '0_0_F_X.'$ct'.bed'
	Rscript expect_vs_obs.R '0_0_F_X.'$ct'.bed' 0_0_F_X.01.bed $ep300_t 2725521370 '0_0_F_X.'$ct'.VAF'
	bedtools window -a $ep300_t -b V_A_0_X.01.bed -w 150 -u > 'V_A_0_X.'$ct'.bed'
	Rscript expect_vs_obs.R 'V_A_0_X.'$ct'.bed' V_A_0_X.01.bed $ep300_t 2725521370 'V_A_0_X.'$ct'.VAF'
	bedtools window -a $ep300_t -b V_0_F_X.01.bed -w 150 -u > 'V_0_F_X.'$ct'.bed'
	Rscript expect_vs_obs.R 'V_0_F_X.'$ct'.bed' V_0_F_X.01.bed $ep300_t 2725521370 'V_0_F_X.'$ct'.VAF'
	bedtools window -a $ep300_t -b 0_A_F_X.01.bed -w 150 -u > '0_A_F_X.'$ct'.bed'
	Rscript expect_vs_obs.R '0_A_F_X.'$ct'.bed' 0_A_F_X.01.bed $ep300_t 2725521370 '0_A_F_X.'$ct'.VAF'
	bedtools window -a $ep300_t -b V_A_F_X.01.bed -w 150 -u > 'V_A_F_X.'$ct'.bed'
	Rscript expect_vs_obs.R 'V_A_F_X.'$ct'.bed' V_A_F_X.01.bed $ep300_t 2725521370 'V_A_F_X.'$ct'.VAF'
	mkdir 'EP300_VAF_'$ct
	mv *'_'*'_'*'_'*'.'$ct'.bed' 'EP300_VAF_'$ct
done

mkdir VAF_enrich
mv *.VAF.enrichment.txt VAF_enrich

echo get heatmap
cd VAF_enrich
ls *VAF.enrichment.txt | awk -F '.' '{print $1}' | sort -u > intersect_pattern.txt
for pattern in $(cat intersect_pattern.txt)
do
	echo $pattern
	tail -n+2 $pattern'.CH12.VAF.enrichment.txt' >> CH12.enrichment.txt
	tail -n+2 $pattern'.Fetal_Liver.VAF.enrichment.txt' >> Fetal_Liver.enrichment.txt
	tail -n+2 $pattern'.MEL.VAF.enrichment.txt' >> MEL.enrichment.txt
done
paste intersect_pattern.txt CH12.enrichment.txt Fetal_Liver.enrichment.txt MEL.enrichment.txt > VAF.enrich.matrix.txt
rm CH12.enrichment.txt Fetal_Liver.enrichment.txt MEL.enrichment.txt

echo plot heatmap
Rscript /Users/universe/Documents/2018_BG/great_analysis_vision_vs_others/windows/get_enrichment_heatmap.R VAF.enrich.matrix.txt VAF.enrich.matrix.png

