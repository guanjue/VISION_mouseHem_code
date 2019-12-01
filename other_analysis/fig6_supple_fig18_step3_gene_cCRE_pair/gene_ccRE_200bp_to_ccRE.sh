#PBS -l nodes=1:ppn=8
#PBS -l walltime=20:00:00
#PBS -j oe
#PBS -A yzz2_e_g_sc_default
#PBS -l pmem=16gb

module load gcc/5.3.1
module load python/2.7.14-anaconda5.0.1
module load bedtools

cd /storage/home/g/gzx103/group/projects/vision/rna/gene_ccRE_200bp_to_ccRE
### sort all mat
#sort -k3,3 vision_rna.gene_ccRE.selected.all.reprod_count120718_withNames.txt > vision_rna.gene_ccRE.selected.all.reprod_count120718_withNames.sort.txt
#sort -k3,3 vision_rna.gene_ccRE.selected.all.reprod_count.WithName1209.txt > vision_rna.gene_ccRE.selected.all.reprod_count.WithName.sort.txt
#sort -k3,3 vision_rna.gene_ccRE.selected.all.reprod_count.WithName1209.txt > vision_rna.gene_ccRE.selected.all.reprod_count.WithName.sort.txt
cp /storage/home/g/gzx103/group/projects/vision/rna/vision_rna.gene_ccRE.selected.all.reprod_count.WithNameSorted.txt vision_rna.gene_ccRE.selected.all.reprod_count.WithName.sort.txt
### extract ccRE-tss pairs
sort -k1,1 -k2,2n vision_cres.bed > vision_cres.sort.bed
time python 200bp_to_ccRE.py

### uniq all info
time sort -u all_converted.txt | sort -k3,3 > all_converted.sort.txt

### get group
time python 200bp_to_ccRE.geneGroup.py

### get cCRE state
cat all_converted.sort.g.txt | awk -F '\t' -v OFS='\t' '{print $1"_"$4"_"$5}' > allcCRE_list.txt
cat /storage/home/gzx103/group/projects/vision/snapshot20_reproduce_16lim_pre_0state_ccRE_is_overestimated/atac_20cell.fun.txt | awk -F '\t' -v OFS='\t' '{print $1"_"$2"_"$3, $5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20,$21,$22,$23,$24}' > atac_fun.matrix.txt
time python vlookup.py -t atac_fun.matrix.txt -m 1 -s allcCRE_list.txt -n 1 -o all_converted.sort.g.cCRE_state.txt -k F
paste all_converted.sort.g.txt all_converted.sort.g.cCRE_state.txt > all_converted.sort.g.cCRE_state.txt.tmp && mv all_converted.sort.g.cCRE_state.txt.tmp all_converted.sort.g.cCRE_state.txt

### split by gene group
cat all_converted.sort.g.txt | awk -F '\t' -v OFS='\t' '{if ($8==1) print $0}' > all_converted.sort.g1.txt
cat all_converted.sort.g.txt | awk -F '\t' -v OFS='\t' '{if ($8==2) print $0}' > all_converted.sort.g2.txt
cat all_converted.sort.g.txt | awk -F '\t' -v OFS='\t' '{if ($8==3) print $0}' > all_converted.sort.g3.txt
cat all_converted.sort.g.txt | awk -F '\t' -v OFS='\t' '{if ($8==4) print $0}' > all_converted.sort.g4.txt

### get uniq
time python 200bp_to_ccRE_uniq.py

### get counts
time python 200bp_to_ccRE_count.py

### plot counts
time Rscript plot_count.R
