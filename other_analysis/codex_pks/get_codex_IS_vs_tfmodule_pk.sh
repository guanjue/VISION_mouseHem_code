##################
### jarkard index
for file in $(cat IS_list.txt)
do
echo $file
i=$(echo $file | awk -F '.' '{print $1}')
echo $i
time sort -k1,1 -k2,2n '/storage/home/gzx103/group/projects/vision/snapshot18_reproduce_0_16lim_all_NB88/index_set_bed/'$file > IS_tmp.sort.bed
### intersect
for j in {1..20}
do
echo $j
time bedtools jaccard -a IS_tmp.sort.bed -b 'tf_tf_network/cluster_'$j'_merged.sort.bed' > 'jaccard_index_'$i'/IS'$i'.tfm'$j'.JI.txt' 
done
done


##################
### Jaccard index matrix
mkdir JI_matrix_folder_IS_vs_TFmodule 
for file in $(cat IS_list.txt)
do
echo $file
i=$(echo $file | awk -F '.' '{print $1}')
echo $i
for j in {1..20}
do
echo $j
time tail -n+2 'jaccard_index_'$i'/IS'$i'.tfm'$j'.JI.txt' >> 'JI_matrix_IS_vs_TFmodule.'$i'.txt'
done
mv 'JI_matrix_IS_vs_TFmodule.'$i'.txt' JI_matrix_folder_IS_vs_TFmodule
done


##################
### get IS index matrix
cat IS_list.txt | awk -F '.' -v OFS='\t' '{print $2}' | awk -F '_' -v OFS='\t' '{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18}' > IS_index_mat.txt
cat IS_list.txt | awk -F '.' -v OFS='\t' '{print $2}' > IS_mat.txt
cp /storage/home/gzx103/scratch/vision/5end/index_set_20cell_pknorm_qda/snapshot18_reproduce_0_16lim_all/atac_20cell.sig.txt IS_sig_mat.txt



LSK	HPC7	CMP	MEP	G1E	ER4	CFUE	ERY	ERY_fl	CFUMK	iMK	GMP	MON	NEU	NK	B	T_CD4	T_CD8
