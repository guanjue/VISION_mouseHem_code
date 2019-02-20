### get ccRE bed files
time cat vision_cres.txt | awk -F ' ' -v OFS='\t' '{print $1,sprintf("%.0f", $2),sprintf("%.0f", $3)}' | sort -k1,1 -k2,2n > vision_cres.bed

### get 200-bp-ccRE bed files
time cat vision_rna.gene_ccRE.selected.all.reprod_count.txt | awk -F '\t' -v OFS='\t' '{print $1,sprintf("%.0f", $5-200),sprintf("%.0f", $5)}' | sort -u | sort -k1,1 -k2,2n > vision_cres_200bp.bed

### get ccRE-to-200bp-bin matched table
time bedtools window -a vision_cres_200bp.bed -b vision_cres.bed -w 1000 > vision_cres_200bp_ccREs_table.bed

cut -f1,2,3 vision_cres_200bp_ccREs_table.bed | sort -u > vision_cres_200bp_ccREs_table.pk.txt

wc -l vision_cres_200bp_ccREs_table.pk.txt


bedtools map -a vision_cres_200bp.bed -b vision_cres.bed -o count > vision_cres_200bp_ccREs_table.1.bed

time sort -u vision_cres_200bp.bed > vision_cres_200bp.uniq.bed

vision_cres_200bp_ccREs_table.bed


time python 200bp_to_ccRE.geneGroup.py

cat all_converted.sort.g.txt | awk -F '\t' -v OFS='\t' '{if ($8==1) print $0}' > all_converted.sort.g1.txt
cat all_converted.sort.g.txt | awk -F '\t' -v OFS='\t' '{if ($8==2) print $0}' > all_converted.sort.g2.txt
cat all_converted.sort.g.txt | awk -F '\t' -v OFS='\t' '{if ($8==3) print $0}' > all_converted.sort.g3.txt
cat all_converted.sort.g.txt | awk -F '\t' -v OFS='\t' '{if ($8==4) print $0}' > all_converted.sort.g4.txt

