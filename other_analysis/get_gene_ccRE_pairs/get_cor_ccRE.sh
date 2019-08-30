cd /storage/home/g/gzx103/group/projects/vision/rna/get_ccRE_gene_pair_nofilter/
### get tss.bed tss 200bp
#head -1 /storage/home/g/gzx103/group/projects/vision/rna/rnaTPM.txt | awk -F ' ' -v OFS='\t' '{print $1, $2, $2+1, $3, $4, $5"_"$6"_"$7"_"$8"_"$9"_"$10"_"$11"_"$12"_"$13"_"$14"_"$15"_"$16}' > gene_tss.bed
tail -n+2 /storage/home/g/gzx103/group/projects/vision/rna/rnaTPM.txt | awk -F ' ' -v OFS='\t' '{print $1, $2, $2+1, $3, $4, $5"_"$6"_"$7"_"$8"_"$9"_"$10"_"$11"_"$12"_"$13"_"$14"_"$15"_"$16}' | sort -u | sort -k1,1 -k2,2n > gene_tss.bed
tail -n+2 /storage/home/g/gzx103/group/projects/vision/pknorm_2_16lim_ref1mo_sample1mo_0424/pknorm_2_16lim_ref1mo_0424_lesshet.state | awk -F ' ' -v OFS='\t' '{print $2,$3,$4, $6"_"$7"_"$9"_"$12"_"$14"_"$18"_"$16"_"$17"_"$20"_"$21"_"$10"_"$13}' | sort -k1,1 -k2,2n > gene_IDEAS_states.bed
### 
time bedtools intersect -a gene_IDEAS_states.bed -b /storage/home/g/gzx103/group/projects/vision/rna/vision_cres.bed -wa > gene_IDEAS_states.inccRE.bed
time bedtools window -a gene_tss.bed -b gene_IDEAS_states.inccRE.bed -w 0 > gene_tss.IDEAS_states.bed
cut -f6 gene_tss.IDEAS_states.bed | awk -F '_' -v OFS='\t' '{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12}' > gene_tss.IDEAS_states.tpm.txt
cut -f10 gene_tss.IDEAS_states.bed | awk -F '_' -v OFS='\t' '{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12}' > gene_tss.IDEAS_states.IDEAS_state.txt

### get tss.bed tss 2MB
tail -n+2 /storage/home/g/gzx103/group/projects/vision/rna/rnaTPM.txt | awk -F ' ' -v OFS='\t' '{if ($2-1000000>0) print $1, $2-1000000, $2+1000000, $3, $4, $5"_"$6"_"$7"_"$8"_"$9"_"$10"_"$11"_"$12"_"$13"_"$14"_"$15"_"$16; else print $1, 0, $2+1000000, $3, $4, $5"_"$6"_"$7"_"$8"_"$9"_"$10"_"$11"_"$12"_"$13"_"$14"_"$15"_"$16}' | sort -u | sort -k1,1 -k2,2n > gene_tss.2MB.bed
### 
time bedtools window -a gene_tss.2MB.bed -b gene_IDEAS_states.inccRE.bed -w 0 > gene_tss.2MB.IDEAS_states.bed
cut -f6 gene_tss.2MB.IDEAS_states.bed | awk -F '_' -v OFS='\t' '{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12}' > gene_tss.2MB.IDEAS_states.tpm.txt
cut -f10 gene_tss.2MB.IDEAS_states.bed | awk -F '_' -v OFS='\t' '{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12}' > gene_tss.IDEAS_states.2MB.IDEAS_state.txt
cut -f1,2,3,4,7,8,9 gene_tss.2MB.IDEAS_states.bed | awk -F '\t' -v OFS='\t' '{print $1"_"$2"_"$3"_"$4, $5"_"$6"_"$7}' > gene_ccRE_pair_table.txt

time shuf -n 1000000 gene_tss.2MB.IDEAS_states.bed > gene_tss.2MB.IDEAS_states.rand.txt
cut -f6 gene_tss.2MB.IDEAS_states.rand.txt | awk -F '_' -v OFS='\t' '{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12}' > gene_tss.2MB.IDEAS_states.tpm.rand.txt
cut -f10 gene_tss.2MB.IDEAS_states.rand.txt | awk -F '_' -v OFS='\t' '{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12}' > gene_tss.IDEAS_states.2MB.IDEAS_state.rand.txt


for i in {1..10}
do
echo $i
### get tss.bed tss 2MB
bin_size=1000
exp_win1=$(($i*$bin_size-1*$bin_size))
exp_win2=$(($i*$bin_size))
#time Rscript expand_bed.R /storage/home/g/gzx103/group/projects/vision/rna/rnaTPM.txt $exp_win 'gene_tss.2MB.'$i'.bed'
time Rscript expand_bed_bin.R /storage/home/g/gzx103/group/projects/vision/rna/rnaTPM.txt $exp_win1 $exp_win2 'gene_tss.2MB.'$i'.bed'
### 
time bedtools window -a 'gene_tss.2MB.'$i'.bed' -b gene_IDEAS_states.inccRE.bed -w 0 > 'gene_tss.2MB.IDEAS_states.'$i'.bed'
cut -f6 'gene_tss.2MB.IDEAS_states.'$i'.bed' | awk -F '_' -v OFS='\t' '{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12}' > 'gene_tss.2MB.IDEAS_states.tpm.'$i'.txt'
cut -f10 'gene_tss.2MB.IDEAS_states.'$i'.bed' | awk -F '_' -v OFS='\t' '{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12}' > 'gene_tss.IDEAS_states.2MB.IDEAS_state.'$i'.txt'
#cut -f1,2,3,4,7,8,9 'gene_tss.2MB.IDEAS_states.'$i'.bed' | awk -F '\t' -v OFS='\t' '{print $1"_"$2"_"$3"_"$4, $5"_"$6"_"$7}' > 'gene_ccRE_pair_table.'$i'.txt'
###
time shuf -n 1000000 'gene_tss.2MB.IDEAS_states.'$i'.bed' > 'gene_tss.2MB.IDEAS_states.rand.'$i'.txt'
cut -f6 'gene_tss.2MB.IDEAS_states.rand.'$i'.txt' | awk -F '_' -v OFS='\t' '{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12}' > 'gene_tss.2MB.IDEAS_states.tpm.rand.'$i'.txt'
cut -f10 'gene_tss.2MB.IDEAS_states.rand.'$i'.txt' | awk -F '_' -v OFS='\t' '{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12}' > 'gene_tss.IDEAS_states.2MB.IDEAS_state.rand.'$i'.txt'
###
done


### get eRP0
time Rscript get_coefficient.R

time Rscript get_nofilter_cor.R




