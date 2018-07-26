echo V_A_F_S
echo V_A_F_S: intersection of VISION and Amit ccREs

echo level uniq
echo V uniq
bedtools window -a atac_20cell.fun.no0.bed.txt -b AmitFriedman_BloodCellEnhancerCatalog_2014_mm10.bed.txt -w 150 -v > V_0_X_X.bed
bedtools window -a V_0_X_X.bed -b FLD_PEAKS_merged_mm10.txt -w 150 -v > V_0_0_X.bed
bedtools window -a V_0_0_X.bed -b SCREEN_FetLiv14_5_mm10_bed_2017.txt -w 150 -v > V_0_0_0.bed
echo A uniq
bedtools window -a AmitFriedman_BloodCellEnhancerCatalog_2014_mm10.bed.txt -b atac_20cell.fun.no0.bed.txt -w 150 -v > 0_A_X_X.bed
bedtools window -a 0_A_X_X.bed -b FLD_PEAKS_merged_mm10.txt -w 150 -v > 0_A_0_X.bed
bedtools window -a 0_A_0_X.bed -b SCREEN_FetLiv14_5_mm10_bed_2017.txt -w 150 -v > 0_A_0_0.bed
echo F uniq
bedtools window -a FLD_PEAKS_merged_mm10.txt -b atac_20cell.fun.no0.bed.txt -w 150 -v > 0_X_F_X.bed
bedtools window -a 0_X_F_X.bed -b AmitFriedman_BloodCellEnhancerCatalog_2014_mm10.bed.txt -w 150 -v > 0_0_F_X.bed
bedtools window -a 0_0_F_X.bed -b SCREEN_FetLiv14_5_mm10_bed_2017.txt -w 150 -v > 0_0_F_0.bed
echo S uniq
bedtools window -a SCREEN_FetLiv14_5_mm10_bed_2017.txt -b atac_20cell.fun.no0.bed.txt -w 150 -v > 0_X_X_S.bed
bedtools window -a 0_X_X_S.bed -b AmitFriedman_BloodCellEnhancerCatalog_2014_mm10.bed.txt -w 150 -v > 0_0_X_S.bed
bedtools window -a 0_0_X_S.bed -b FLD_PEAKS_merged_mm10.txt -w 150 -v > 0_0_0_S.bed


echo level 2 uniq
echo V_A_0_0 
bedtools window -a atac_20cell.fun.no0.bed.txt -b AmitFriedman_BloodCellEnhancerCatalog_2014_mm10.bed.txt -w 150 -u > V_A_X_X.bed
bedtools window -a V_A_X_X.bed -b FLD_PEAKS_merged_mm10.txt -w 150 -v > V_A_0_X.bed
bedtools window -a V_A_0_X.bed -b SCREEN_FetLiv14_5_mm10_bed_2017.txt -w 150 -v > V_A_0_0.bed
echo V_0_F_0 
bedtools window -a V_0_X_X.bed -b FLD_PEAKS_merged_mm10.txt -w 150 -u > V_0_F_X.bed
bedtools window -a V_0_F_X.bed -b SCREEN_FetLiv14_5_mm10_bed_2017.txt -w 150 -v > V_0_F_0.bed
echo V_0_0_S 
bedtools window -a V_0_0_X.bed -b SCREEN_FetLiv14_5_mm10_bed_2017.txt -w 150 -u > V_0_0_S.bed
echo 0_A_F_0 
bedtools window -a AmitFriedman_BloodCellEnhancerCatalog_2014_mm10.bed.txt -b FLD_PEAKS_merged_mm10.txt -w 150 -u > X_A_F_X.bed
bedtools window -a X_A_F_X.bed -b atac_20cell.fun.no0.bed.txt -w 150 -v > 0_A_F_X.bed
bedtools window -a 0_A_F_X.bed -b SCREEN_FetLiv14_5_mm10_bed_2017.txt -w 150 -v > 0_A_F_0.bed
echo 0_A_0_S 
bedtools window -a AmitFriedman_BloodCellEnhancerCatalog_2014_mm10.bed.txt -b SCREEN_FetLiv14_5_mm10_bed_2017.txt -w 150 -u > X_A_X_S.bed
bedtools window -a X_A_X_S.bed -b atac_20cell.fun.no0.bed.txt -w 150 -v > 0_A_X_S.bed
bedtools window -a 0_A_X_S.bed -b FLD_PEAKS_merged_mm10.txt -w 150 -v > 0_A_0_S.bed
echo 0_0_F_S 
bedtools window -a FLD_PEAKS_merged_mm10.txt -b SCREEN_FetLiv14_5_mm10_bed_2017.txt -w 150 -u > X_X_F_S.bed
bedtools window -a X_X_F_S.bed -b atac_20cell.fun.no0.bed.txt -w 150 -v > 0_X_F_S.bed
bedtools window -a 0_X_F_S.bed -b AmitFriedman_BloodCellEnhancerCatalog_2014_mm10.bed.txt -w 150 -v > 0_0_F_S.bed

echo level 3 uniq
echo V_A_F_0 
bedtools window -a V_A_X_X.bed -b FLD_PEAKS_merged_mm10.txt -w 150 -u > V_A_F_X.bed
bedtools window -a V_A_F_X.bed -b SCREEN_FetLiv14_5_mm10_bed_2017.txt -w 150 -v > V_A_F_0.bed
echo V_A_0_S 
bedtools window -a V_A_X_X.bed -b SCREEN_FetLiv14_5_mm10_bed_2017.txt -w 150 -u > V_A_X_S.bed
bedtools window -a V_A_X_S.bed -b FLD_PEAKS_merged_mm10.txt -w 150 -v > V_A_0_S.bed
echo V_0_F_S 
bedtools window -a atac_20cell.fun.no0.bed.txt -b FLD_PEAKS_merged_mm10.txt -w 150 -u > V_X_F_X.bed
bedtools window -a V_X_F_X.bed -b SCREEN_FetLiv14_5_mm10_bed_2017.txt -w 150 -u > V_X_F_S.bed
bedtools window -a V_X_F_S.bed -b AmitFriedman_BloodCellEnhancerCatalog_2014_mm10.bed.txt -w 150 -v > V_0_F_S.bed
echo 0_A_F_S
bedtools window -a X_A_F_X.bed -b SCREEN_FetLiv14_5_mm10_bed_2017.txt -w 150 -u > X_A_F_S.bed
bedtools window -a X_A_F_S.bed -b atac_20cell.fun.no0.bed.txt -w 150 -v > 0_A_F_S.bed

echo level 4 uniq
bedtools window -a V_A_X_X.bed -b FLD_PEAKS_merged_mm10.txt -w 150 -u > V_A_F_X.bed
bedtools window -a V_A_F_X.bed -b SCREEN_FetLiv14_5_mm10_bed_2017.txt -w 150 -u > V_A_F_S.bed








cut -f1,2,3 V_A_F_S.bed > V_A_F_S.01.bed
cut -f1,2,3 V_A_F_0.bed > V_A_F_0.01.bed
cut -f1,2,3 V_A_0_S.bed > V_A_0_S.01.bed
cut -f1,2,3 V_0_F_S.bed > V_0_F_S.01.bed
cut -f1,2,3 0_A_F_S.bed > 0_A_F_S.01.bed
cut -f1,2,3 V_A_0_0.bed > V_A_0_0.01.bed
cut -f1,2,3 V_0_F_0.bed > V_0_F_0.01.bed
cut -f1,2,3 0_A_F_0.bed > 0_A_F_0.01.bed
cut -f1,2,3 V_0_0_S.bed > V_0_0_S.01.bed
cut -f1,2,3 0_A_0_S.bed > 0_A_0_S.01.bed
cut -f1,2,3 0_0_F_S.bed > 0_0_F_S.01.bed
cut -f1,2,3 0_0_0_S.bed > 0_0_0_S.01.bed
cut -f1,2,3 0_0_F_0.bed > 0_0_F_0.01.bed
cut -f1,2,3 0_A_0_0.bed > 0_A_0_0.01.bed
cut -f1,2,3 V_0_0_0.bed > V_0_0_0.01.bed


bedtools window -a Known_erythroid_REs_2017_mm10.txt -b V_A_F_S.01.bed -w 150 -u > V_A_F_S.01.K.bed
bedtools window -a Known_erythroid_REs_2017_mm10.txt -b V_A_F_0.01.bed -w 150 -u > V_A_F_0.01.K.bed
bedtools window -a Known_erythroid_REs_2017_mm10.txt -b V_A_0_S.01.bed -w 150 -u > V_A_0_S.01.K.bed
bedtools window -a Known_erythroid_REs_2017_mm10.txt -b V_0_F_S.01.bed -w 150 -u > V_0_F_S.01.K.bed
bedtools window -a Known_erythroid_REs_2017_mm10.txt -b 0_A_F_S.01.bed -w 150 -u > 0_A_F_S.01.K.bed
bedtools window -a Known_erythroid_REs_2017_mm10.txt -b V_A_0_0.01.bed -w 150 -u > V_A_0_0.01.K.bed
bedtools window -a Known_erythroid_REs_2017_mm10.txt -b V_0_F_0.01.bed -w 150 -u > V_0_F_0.01.K.bed
bedtools window -a Known_erythroid_REs_2017_mm10.txt -b 0_A_F_0.01.bed -w 150 -u > 0_A_F_0.01.K.bed
bedtools window -a Known_erythroid_REs_2017_mm10.txt -b V_0_0_S.01.bed -w 150 -u > V_0_0_S.01.K.bed
bedtools window -a Known_erythroid_REs_2017_mm10.txt -b 0_A_0_S.01.bed -w 150 -u > 0_A_0_S.01.K.bed
bedtools window -a Known_erythroid_REs_2017_mm10.txt -b 0_0_F_S.01.bed -w 150 -u > 0_0_F_S.01.K.bed
bedtools window -a Known_erythroid_REs_2017_mm10.txt -b 0_0_0_S.01.bed -w 150 -u > 0_0_0_S.01.K.bed
bedtools window -a Known_erythroid_REs_2017_mm10.txt -b 0_0_F_0.01.bed -w 150 -u > 0_0_F_0.01.K.bed
bedtools window -a Known_erythroid_REs_2017_mm10.txt -b 0_A_0_0.01.bed -w 150 -u > 0_A_0_0.01.K.bed
bedtools window -a Known_erythroid_REs_2017_mm10.txt -b V_0_0_0.01.bed -w 150 -u > V_0_0_0.01.K.bed






















