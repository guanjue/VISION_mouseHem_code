cd /storage/home/gzx103/group/projects/vision/pknorm_2_16lim_ref1mo_sample1mo_0424

### get IDEAS state bed
time tail -n+2 ~/group/projects/vision/pknorm_2_16lim_ref1mo_sample1mo_0424/pknorm_2_16lim_ref1mo_0424_lesshet.state | awk -F ' ' -v OFS='\t' '{print $2,$3,$4,$5}' | sort -k1,1 -k2,2n > B_SPL.state.bed
time tail -n+2 ~/group/projects/vision/pknorm_2_16lim_ref1mo_sample1mo_0424/pknorm_2_16lim_ref1mo_0424_lesshet.state | awk -F ' ' -v OFS='\t' '{print $2,$3,$4,$6}' | sort -k1,1 -k2,2n > CFU_E_ad.state.bed
time tail -n+2 ~/group/projects/vision/pknorm_2_16lim_ref1mo_sample1mo_0424/pknorm_2_16lim_ref1mo_0424_lesshet.state | awk -F ' ' -v OFS='\t' '{print $2,$3,$4,$7}' | sort -k1,1 -k2,2n > CFUMK.state.bed
time tail -n+2 ~/group/projects/vision/pknorm_2_16lim_ref1mo_sample1mo_0424/pknorm_2_16lim_ref1mo_0424_lesshet.state | awk -F ' ' -v OFS='\t' '{print $2,$3,$4,$8}' | sort -k1,1 -k2,2n > CLP.state.bed
time tail -n+2 ~/group/projects/vision/pknorm_2_16lim_ref1mo_sample1mo_0424/pknorm_2_16lim_ref1mo_0424_lesshet.state | awk -F ' ' -v OFS='\t' '{print $2,$3,$4,$9}' | sort -k1,1 -k2,2n > CMP.state.bed
time tail -n+2 ~/group/projects/vision/pknorm_2_16lim_ref1mo_sample1mo_0424/pknorm_2_16lim_ref1mo_0424_lesshet.state | awk -F ' ' -v OFS='\t' '{print $2,$3,$4,$10}' | sort -k1,1 -k2,2n > ER4.state.bed
time tail -n+2 ~/group/projects/vision/pknorm_2_16lim_ref1mo_sample1mo_0424/pknorm_2_16lim_ref1mo_0424_lesshet.state | awk -F ' ' -v OFS='\t' '{print $2,$3,$4,$11}' | sort -k1,1 -k2,2n > ERY_ad.state.bed
time tail -n+2 ~/group/projects/vision/pknorm_2_16lim_ref1mo_sample1mo_0424/pknorm_2_16lim_ref1mo_0424_lesshet.state | awk -F ' ' -v OFS='\t' '{print $2,$3,$4,$12}' | sort -k1,1 -k2,2n > ERY_fl.state.bed
time tail -n+2 ~/group/projects/vision/pknorm_2_16lim_ref1mo_sample1mo_0424/pknorm_2_16lim_ref1mo_0424_lesshet.state | awk -F ' ' -v OFS='\t' '{print $2,$3,$4,$13}' | sort -k1,1 -k2,2n > G1E.state.bed
time tail -n+2 ~/group/projects/vision/pknorm_2_16lim_ref1mo_sample1mo_0424/pknorm_2_16lim_ref1mo_0424_lesshet.state | awk -F ' ' -v OFS='\t' '{print $2,$3,$4,$14}' | sort -k1,1 -k2,2n > GMP.state.bed
time tail -n+2 ~/group/projects/vision/pknorm_2_16lim_ref1mo_sample1mo_0424/pknorm_2_16lim_ref1mo_0424_lesshet.state | awk -F ' ' -v OFS='\t' '{print $2,$3,$4,$15}' | sort -k1,1 -k2,2n > HPC7.state.bed
time tail -n+2 ~/group/projects/vision/pknorm_2_16lim_ref1mo_sample1mo_0424/pknorm_2_16lim_ref1mo_0424_lesshet.state | awk -F ' ' -v OFS='\t' '{print $2,$3,$4,$16}' | sort -k1,1 -k2,2n > LSK_BM.state.bed
time tail -n+2 ~/group/projects/vision/pknorm_2_16lim_ref1mo_sample1mo_0424/pknorm_2_16lim_ref1mo_0424_lesshet.state | awk -F ' ' -v OFS='\t' '{print $2,$3,$4,$17}' | sort -k1,1 -k2,2n > MEP.state.bed
time tail -n+2 ~/group/projects/vision/pknorm_2_16lim_ref1mo_sample1mo_0424/pknorm_2_16lim_ref1mo_0424_lesshet.state | awk -F ' ' -v OFS='\t' '{print $2,$3,$4,$18}' | sort -k1,1 -k2,2n > MK_imm_ad.state.bed
time tail -n+2 ~/group/projects/vision/pknorm_2_16lim_ref1mo_sample1mo_0424/pknorm_2_16lim_ref1mo_0424_lesshet.state | awk -F ' ' -v OFS='\t' '{print $2,$3,$4,$19}' | sort -k1,1 -k2,2n > MK_mat_fl.state.bed
time tail -n+2 ~/group/projects/vision/pknorm_2_16lim_ref1mo_sample1mo_0424/pknorm_2_16lim_ref1mo_0424_lesshet.state | awk -F ' ' -v OFS='\t' '{print $2,$3,$4,$20}' | sort -k1,1 -k2,2n > MONO_BM.state.bed
time tail -n+2 ~/group/projects/vision/pknorm_2_16lim_ref1mo_sample1mo_0424/pknorm_2_16lim_ref1mo_0424_lesshet.state | awk -F ' ' -v OFS='\t' '{print $2,$3,$4,$21}' | sort -k1,1 -k2,2n > NEU.state.bed
time tail -n+2 ~/group/projects/vision/pknorm_2_16lim_ref1mo_sample1mo_0424/pknorm_2_16lim_ref1mo_0424_lesshet.state | awk -F ' ' -v OFS='\t' '{print $2,$3,$4,$22}' | sort -k1,1 -k2,2n > NK_SPL.state.bed
time tail -n+2 ~/group/projects/vision/pknorm_2_16lim_ref1mo_sample1mo_0424/pknorm_2_16lim_ref1mo_0424_lesshet.state | awk -F ' ' -v OFS='\t' '{print $2,$3,$4,$23}' | sort -k1,1 -k2,2n > T_CD4_SPL.state.bed
time tail -n+2 ~/group/projects/vision/pknorm_2_16lim_ref1mo_sample1mo_0424/pknorm_2_16lim_ref1mo_0424_lesshet.state | awk -F ' ' -v OFS='\t' '{print $2,$3,$4,$24}' | sort -k1,1 -k2,2n > T_CD8_SPL.state.bed

### get ccRE list
cat ~/group/projects/vision/rna/vision_cres.txt | awk -F ' ' -v OFS='\t' '{print $1, $2, $3, $1"_"$2"_"$3}' > vision_cres.215117.bed

### get split bed intersect
for ct in $(cat ct_list.txt)
do
time bedtools intersect -a vision_cres.215117.bed -b $ct'.state.bed' > 'ccRE.'$ct'.state.1.bed'
time sort -k1,1 -k2,2n 'ccRE.'$ct'.state.1.bed' > 'ccRE.'$ct'.state.1.sort.bed'
time bedtools map -a 'ccRE.'$ct'.state.1.sort.bed' -b $ct'.state.bed' -c 4 -o collapse > 'ccRE.'$ct'.state.2.bed'
time sort -k4,4 'ccRE.'$ct'.state.2.bed' > 'ccRE.'$ct'.state.2.pksort.bed'
#rm $ct'.state.bed'
rm 'ccRE.'$ct'.state.2.bed'
rm 'ccRE.'$ct'.state.1.bed'
rm 'ccRE.'$ct'.state.1.sort.bed'
done


for ct in $(cat ct_list.txt)
do
### split for speed up
mkdir $ct'_folder'
for i in {0..61}
do
	echo $i
	rownum=$((i * 10000 + 1))
	tail -n+$rownum 'ccRE.'$ct'.state.2.pksort.bed' | head -10000 > 'ccRE.'$ct'.state.2.'$i'.bed'
	mv 'ccRE.'$ct'.state.2.'$i'.bed' $ct'_folder'
done
### get IDEAS state percentage tracks
for i in {0..61}
do
echo \#PBS -l nodes=1:ppn=8 > 'PBS.'$i'.sh'
echo \#PBS -l walltime=10:00:00 >> 'PBS.'$i'.sh'
echo \#PBS -j oe >> 'PBS.'$i'.sh'
echo \#PBS -A yzz2_e_g_sc_default >> 'PBS.'$i'.sh'
echo \#PBS -l pmem=16gb >> 'PBS.'$i'.sh'
echo module load gcc/5.3.1 >> 'PBS.'$i'.sh'
echo module load python/2.7.14-anaconda5.0.1 >> 'PBS.'$i'.sh'
echo module load bedtools >> 'PBS.'$i'.sh'
echo cd /storage/home/gzx103/group/projects/vision/pknorm_2_16lim_ref1mo_sample1mo_0424/get_eRP_tracks/$ct'_folder' >> 'PBS.'$i'.sh' 
echo time Rscript ../get_8tracks.R 'ccRE.ERY_ad.state.2.'$i'.bed' 'ccRE.ERY_ad.state.2.'$i'.statepercent.bed' >> 'PBS.'$i'.sh'
qsub 'PBS.'$i'.sh'
done
done


for ct in $(cat ct_list.txt)
do
cd /storage/home/gzx103/group/projects/vision/pknorm_2_16lim_ref1mo_sample1mo_0424/get_eRP_tracks/$ct'_folder'
### merge all percentage
rm 'ccRE.'$ct'.state.2.all.statepercent.txt'
for i in {0..61}
do
cat 'ccRE.'$ct'.state.2.'$i'.statepercent.bed' >> 'ccRE.'$ct'.state.2.all.statepercent.txt'
done
### merge all peaks
rm 'ccRE.'$ct'.state.2.all.bed'
for i in {0..61}
do
cat 'ccRE.'$ct'.state.2.'$i'.bed' | cut -f4 | uniq | awk -F '_' -v OFS='\t' '{print $1,$2,$3, $1"_"$2"_"$3}' >> 'ccRE.'$ct'.state.2.all.bed'
done
done



for ct in $(cat ct_list.txt)
do
###
cd $ct'_folder'
### get percentage
time Rscript ../get_eRPs.8track.R $ct
### get bigwig
for i in {1..4}
do
	echo $i
	sort -k1,1 -k2,2n $ct'.tss.'$i'.bedgraph' | awk -F '\t' -v OFS='\t' '{if ($4!="NA") print $0}' > $ct'.tss.'$i'.sort.bedgraph'
	bedtools merge -i $ct'.tss.'$i'.sort.bedgraph' -c 4 -o mean > $ct'.tss.'$i'.sort.merge.bedgraph'
	time ~/group/software/ucsc/bedGraphToBigWig $ct'.tss.'$i'.sort.merge.bedgraph' ~/group/genome/mm10/mm10.chrom.sort.sizes $ct'.tss.'$i'.sort.bw'
done
for i in {1..4}
do
	echo $i
	sort -k1,1 -k2,2n $ct'.dist.'$i'.bedgraph' | awk -F '\t' -v OFS='\t' '{if ($4!="NA") print $0}' > $ct'.dist.'$i'.sort.bedgraph'
	bedtools merge -i $ct'.dist.'$i'.sort.bedgraph' -c 4 -o mean > $ct'.dist.'$i'.sort.merge.bedgraph'
	time ~/group/software/ucsc/bedGraphToBigWig $ct'.dist.'$i'.sort.merge.bedgraph' ~/group/genome/mm10/mm10.chrom.sort.sizes $ct'.dist.'$i'.sort.bw'
done
cd ..
done





