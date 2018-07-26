### VandA: intersection of VISION and Amit ccREs
echo VandA: intersection of VISION and Amit ccREs
bedtools intersect -a atac_20cell.fun.no0.bed.txt -b AmitFriedman_BloodCellEnhancerCatalog_2014_mm10.bed.txt -wa > VandA.bed

### VandFandSnotA: ccREs that overlap among VIS, FLD, and SCR but excluding Amit ccREs
echo VandFandSnotA: ccREs that overlap among VIS, FLD, and SCR but excluding Amit ccREs
bedtools intersect -a atac_20cell.fun.no0.bed.txt -b SCREEN_FetLiv14_5_mm10_bed_2017.txt -wa > VandS.bed
bedtools intersect -a VandS.bed -b FLD_Merge_DNase_peaks_for_all_4_cell_types_E14.5_FL_StamBender.bed.txt -wa > VandSF.bed
bedtools intersect -a VandSF.bed -b AmitFriedman_BloodCellEnhancerCatalog_2014_mm10.bed.txt -v > VandSFnotA.bed
rm VandS.bed
rm VandSF.bed

### VIS_ccREs_not_SFA: VISION ccREs that are not in any of SCR, FLD or Amit
echo VIS_ccREs_not_SFA: VISION ccREs that are not in any of SCR, FLD or Amit
bedtools intersect -a atac_20cell.fun.no0.bed.txt -b SCREEN_FetLiv14_5_mm10_bed_2017.txt -v > VnotS.bed
bedtools intersect -a VnotS.bed -b FLD_Merge_DNase_peaks_for_all_4_cell_types_E14.5_FL_StamBender.bed.txt -v > VnotSF.bed
bedtools intersect -a VnotSF.bed -b AmitFriedman_BloodCellEnhancerCatalog_2014_mm10.bed.txt -v > VnotSFA.bed
rm VnotS.bed
rm VnotSF.bed

### FLD_not_SVA: Fetal liver DNaseI HSs not in SCR, VIS, or Amit
echo FLD_not_SVA: Fetal liver DNaseI HSs not in SCR, VIS, or Amit
bedtools intersect -a FLD_Merge_DNase_peaks_for_all_4_cell_types_E14.5_FL_StamBender.bed.txt -b SCREEN_FetLiv14_5_mm10_bed_2017.txt -v > FLDnotS.bed
bedtools intersect -a FLDnotS.bed -b atac_20cell.fun.no0.bed.txt -v > FLDnotSV.bed
bedtools intersect -a FLDnotSV.bed -b AmitFriedman_BloodCellEnhancerCatalog_2014_mm10.bed.txt -v > FLDnotSvA.bed
rm FLDnotS.bed
rm FLDnotSV.bed

### SCRnotVFA: SCREEN ccREs not in VIS, FLD, or Amit
echo SCRnotVFA: SCREEN ccREs not in VIS, FLD, or Amit
bedtools intersect -a SCREEN_FetLiv14_5_mm10_bed_2017.txt -b atac_20cell.fun.no0.bed.txt -v > SnotV.bed
bedtools intersect -a SnotV.bed -b FLD_Merge_DNase_peaks_for_all_4_cell_types_E14.5_FL_StamBender.bed.txt -v > SnotVF.bed
bedtools intersect -a SnotVF.bed -b AmitFriedman_BloodCellEnhancerCatalog_2014_mm10.bed.txt -v > SnotVFA.bed
rm SnotV.bed
rm SnotVF.bed


#atac_20cell.fun.no0.bed.txt
#SCREEN_FetLiv14_5_mm10_bed_2017.txt
#FLD_Merge_DNase_peaks_for_all_4_cell_types_E14.5_FL_StamBender.bed.txt
#AmitFriedman_BloodCellEnhancerCatalog_2014_mm10.bed.txt



cut -f1,2,3 FLDnotSvA.bed > FLDnotSvA_01.bed
cut -f1,2,3 SnotVFA.bed > SnotVFA_01.bed
cut -f1,2,3 VandA.bed > VandA_01.bed
cut -f1,2,3 VandSFnotA.bed > VandSFnotA_01.bed
cut -f1,2,3 VnotSFA.bed > VnotSFA_01.bed






bedtools intersect -a Known_erythroid_REs_2017_mm10.txt -b FLDnotSvA_01.bed -wa > FLDnotSvA_01_K.bed
bedtools intersect -a Known_erythroid_REs_2017_mm10.txt -b SnotVFA_01.bed -wa > SnotVFA_01_K.bed
bedtools intersect -a Known_erythroid_REs_2017_mm10.txt -b VandA_01.bed -wa > VandA_01_K.bed
bedtools intersect -a Known_erythroid_REs_2017_mm10.txt -b VandSFnotA_01.bed -wa > VandSFnotA_01_K.bed
bedtools intersect -a Known_erythroid_REs_2017_mm10.txt -b VnotSFA_01.bed -wa > VnotSFA_01_K.bed






















