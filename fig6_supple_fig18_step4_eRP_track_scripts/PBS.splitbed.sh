#PBS -l nodes=1:ppn=8
#PBS -l walltime=10:00:00
#PBS -j oe
#PBS -A yzz2_e_g_sc_default
#PBS -l pmem=16gb
module load gcc/5.3.1
module load python/2.7.14-anaconda5.0.1
module load bedtools
cd /storage/home/gzx103/group/projects/vision/pknorm_2_16lim_ref1mo_sample1mo_0424/get_eRP_tracks

state_file=/storage/home/gzx103/group/projects/vision/pknorm_2_16lim_ref1mo_sample1mo_0424/pknorm_2_16lim_ref1mo_0424_lesshet.state
ct_list_file=/storage/home/gzx103/group/projects/vision/pknorm_2_16lim_ref1mo_sample1mo_0424/get_eRP_tracks/ct_list.txt

### get state of each file
time tail -n+2 $state_file | awk -F ' ' -v OFS='\t' '{print $2,$3,$4,$5}' | sort -k1,1 -k2,2n > B_SPL.state.bed
time tail -n+2 $state_file | awk -F ' ' -v OFS='\t' '{print $2,$3,$4,$6}' | sort -k1,1 -k2,2n > CFU_E_ad.state.bed
time tail -n+2 $state_file | awk -F ' ' -v OFS='\t' '{print $2,$3,$4,$7}' | sort -k1,1 -k2,2n > CFUMK.state.bed
time tail -n+2 $state_file | awk -F ' ' -v OFS='\t' '{print $2,$3,$4,$8}' | sort -k1,1 -k2,2n > CLP.state.bed
time tail -n+2 $state_file | awk -F ' ' -v OFS='\t' '{print $2,$3,$4,$9}' | sort -k1,1 -k2,2n > CMP.state.bed
time tail -n+2 $state_file | awk -F ' ' -v OFS='\t' '{print $2,$3,$4,$10}' | sort -k1,1 -k2,2n > ER4.state.bed
time tail -n+2 $state_file | awk -F ' ' -v OFS='\t' '{print $2,$3,$4,$11}' | sort -k1,1 -k2,2n > ERY_ad.state.bed
time tail -n+2 $state_file | awk -F ' ' -v OFS='\t' '{print $2,$3,$4,$12}' | sort -k1,1 -k2,2n > ERY_fl.state.bed
time tail -n+2 $state_file | awk -F ' ' -v OFS='\t' '{print $2,$3,$4,$13}' | sort -k1,1 -k2,2n > G1E.state.bed
time tail -n+2 $state_file | awk -F ' ' -v OFS='\t' '{print $2,$3,$4,$14}' | sort -k1,1 -k2,2n > GMP.state.bed
time tail -n+2 $state_file | awk -F ' ' -v OFS='\t' '{print $2,$3,$4,$15}' | sort -k1,1 -k2,2n > HPC7.state.bed
time tail -n+2 $state_file | awk -F ' ' -v OFS='\t' '{print $2,$3,$4,$16}' | sort -k1,1 -k2,2n > LSK_BM.state.bed
time tail -n+2 $state_file | awk -F ' ' -v OFS='\t' '{print $2,$3,$4,$17}' | sort -k1,1 -k2,2n > MEP.state.bed
time tail -n+2 $state_file | awk -F ' ' -v OFS='\t' '{print $2,$3,$4,$18}' | sort -k1,1 -k2,2n > MK_imm_ad.state.bed
time tail -n+2 $state_file | awk -F ' ' -v OFS='\t' '{print $2,$3,$4,$19}' | sort -k1,1 -k2,2n > MK_mat_fl.state.bed
time tail -n+2 $state_file | awk -F ' ' -v OFS='\t' '{print $2,$3,$4,$20}' | sort -k1,1 -k2,2n > MONO_BM.state.bed
time tail -n+2 $state_file | awk -F ' ' -v OFS='\t' '{print $2,$3,$4,$21}' | sort -k1,1 -k2,2n > NEU.state.bed
time tail -n+2 $state_file | awk -F ' ' -v OFS='\t' '{print $2,$3,$4,$22}' | sort -k1,1 -k2,2n > NK_SPL.state.bed
time tail -n+2 $state_file | awk -F ' ' -v OFS='\t' '{print $2,$3,$4,$23}' | sort -k1,1 -k2,2n > T_CD4_SPL.state.bed
time tail -n+2 $state_file | awk -F ' ' -v OFS='\t' '{print $2,$3,$4,$24}' | sort -k1,1 -k2,2n > T_CD8_SPL.state.bed

### get split bed intersect
for ct in $(cat $ct_list_file)
do
	time bedtools intersect -a vision_cres.215117.bed -b $ct'.state.bed' > 'ccRE.'$ct'.state.1.bed'
	#cat atac_20cell.sort.bed | awk -F '\t' -v OFS='\t' '{print $3-$2}' > ccRE.length.txt
	time sort -k1,1 -k2,2n 'ccRE.'$ct'.state.1.bed' > 'ccRE.'$ct'.state.1.sort.bed'
	time bedtools map -a 'ccRE.'$ct'.state.1.sort.bed' -b $ct'.state.bed' -c 4 -o collapse > 'ccRE.'$ct'.state.2.bed'
	time sort -k4,4 'ccRE.'$ct'.state.2.bed' > 'ccRE.'$ct'.state.2.pksort.bed'
	#rm $ct'.state.bed'
	rm 'ccRE.'$ct'.state.2.bed'
	rm 'ccRE.'$ct'.state.1.bed'
	rm 'ccRE.'$ct'.state.1.sort.bed'
done


