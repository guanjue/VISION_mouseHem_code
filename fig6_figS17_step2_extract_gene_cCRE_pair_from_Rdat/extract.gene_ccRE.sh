#PBS -l nodes=1:ppn=8
#PBS -l walltime=20:00:00
#PBS -j oe
#PBS -A yzz2_e_g_sc_default
#PBS -l pmem=16gb

module load gcc/5.3.1
module load python/2.7.14-anaconda5.0.1
module load bedtools

cd /storage/home/gzx103/group/projects/vision/rna

time Rscript get_gene_groups.R

### extract information
for i in {1..12}
do
	echo $i
	cat 'vision_rna_tss2k_ccreunit.chr'*'.'*'.'$i'.gene_ccRE.txt' > 'vision_rna.gene_ccRE.c'$i'.selected.all.txt'
done

### get reproducibility in leave-one-out
time python check_reproducibility.py

###
#time cat gencode.vM4.annotation.gtf | awk -F '"' -v OFS='\t' '{print $1,$2,$10}' | awk -F '\t' -v OFS='\t' '{if ($7=="+") print $1"_"$4, $11":"$10":"$7; else if ($7=="-") print $1"_"$5, $11":"$10":"$7}' > gencode.vM4.tss.tss2gene.txt
#time cat gencode.vM4.annotation.gtf | awk -F '"' -v OFS='\t' '{print $1,$2,$10}' | awk -F '\t' -v OFS='\t' '{if ($7=="+") print $1"_"$4, $11":"$10; else if ($7=="-") print $1"_"$5, $11":"$10}' > gencode.vM4.tss.tss2gene0.txt
###

### get gene name
time python get_gene_nameid.py

### get gene group tss list
for i in {1..4}
do
	echo $i
	cat 'vision_rna_tss2k_ccreunit.chr'*'.'$i'.'*'.gene_ccRE.txt' | awk -F ' ' -v OFS='\t' '{print $1"_"$2}' | sort -u > 'gene.group.'$i'.txt'
done


sort -k3,3 vision_rna.gene_ccRE.selected.all.reprod_count.WithName1209.txt > vision_rna.gene_ccRE.selected.all.reprod_count.WithNameSorted.txt

cp vision_rna.gene_ccRE.selected.all.reprod_count.WithName1209.txt gene_ccRE_200bp_to_ccRE/
cp gene.group.*.txt gene_ccRE_200bp_to_ccRE/


