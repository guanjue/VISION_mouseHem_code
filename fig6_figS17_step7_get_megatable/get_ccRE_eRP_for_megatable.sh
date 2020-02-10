#PBS -l nodes=1:ppn=8
#PBS -l walltime=10:00:00
#PBS -j oe
#PBS -A yzz2_e_g_sc_default
#PBS -l pmem=16gb

module load gcc/5.3.1
module load python/2.7.14-anaconda5.0.1
module load bedtools


cd ~/group/projects/vision/rna/ccRE_eRPs_Rdat_version/

#cp /storage/home/gzx103/group/projects/vision/rna/mouse_allChr_wTAD_022519.txt /storage/home/gzx103/group/projects/vision/rna/ccRE_eRPs_Rdat_version/
#cd /storage/home/gzx103/group/projects/vision/rna/ccRE_eRPs_Rdat_version
#cut -f1,4,5 mouse_allChr_wTAD_022519.txt | awk -F '\t' -v OFS='\t' '{print $1,$2,$3,$1"_"$2"_"$3}' | sort -k1,1 -k2,2n | uniq | tail -n+2 > mouse_allChr_wTAD_022519.bed
#cut -f1,4,5 mouse_allChr_wTAD_022519.txt | awk -F '\t' -v OFS='\t' '{if ($1!="chr") print $1,$2,$3,$1"_"$2"_"$3}' > mouse_allChr_wTAD_022519.repeat.bed
cut -f1,4,5 /storage/home/gzx103/group/projects/vision/rna/gene_ccRE_200bp_to_ccRE/all_converted.sort.g.cCRE_state.txt | awk -F '\t' -v OFS='\t' '{if ($1!="chr") print $1,$2,$3,$1"_"$2"_"$3}' > mouse_allChr_wTAD_022519.repeat.bed

### get ccRE table
#cp ~/group/projects/vision/pknorm_2_16lim_ref1mo_sample1mo_0424/get_eRP_tracks/ccRE.ERY_ad.state.2.all.bed ./

declare -a ct_vector=("B_SPL" "CFU_E_ad" "CFUMK" "CLP" "CMP" "ER4" "ERY_ad" "ERY_fl" "G1E" "GMP" "HPC7" "LSK_BM" "MEP" "MK_imm_ad" "MK_mat_fl" "MONO_BM" "NEU" "NK_SPL" "T_CD4_SPL" "T_CD8_SPL")
declare -a pos_vector=("tss" "dist")

cat /storage/home/gzx103/group/projects/vision/pknorm_2_16lim_ref1mo_sample1mo_0424/get_eRP_tracks/T_CD8_SPL_folder/T_CD8_SPL.dist.1.bedgraph | awk -F '\t' -v OFS='\t' '{print $1"_"$2"_"$3}' > ccRE_eRP.table.txt

for ct in "${ct_vector[@]}"
do
	for pos in "${pos_vector[@]}"
	do
		for t in {1..4}
		do
			cut -f4 /storage/home/gzx103/group/projects/vision/pknorm_2_16lim_ref1mo_sample1mo_0424/get_eRP_tracks/$ct'_folder/'$ct'.'$pos'.'$t'.bedgraph' > tmp.eRP.txt
			paste ccRE_eRP.table.txt tmp.eRP.txt > ccRE_eRP.table.txt.tmp && mv ccRE_eRP.table.txt.tmp ccRE_eRP.table.txt
		done
	done
done

### round to reduce file size
time Rscript round.R
#paste ccRE.ERY_ad.state.2.all.bed tss.1.bedgraph tss.2.bedgraph tss.3.bedgraph tss.4.bedgraph dist.1.bedgraph dist.2.bedgraph dist.3.bedgraph dist.4.bedgraph | cut -f4,8,12,16,20,24,28,32,36 > ccRE_eRP.table.txt


time python get_megatable_eRP.py

cat header.txt mouse_allChr_wTAD_022519.eRPs.txt > mouse_allChr_wTAD_022519.eRPs.txt.tmp && mv mouse_allChr_wTAD_022519.eRPs.txt.tmp mouse_allChr_wTAD_022519.eRPs.txt
