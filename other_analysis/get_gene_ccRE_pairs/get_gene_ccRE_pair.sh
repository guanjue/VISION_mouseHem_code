for i in {1..12}
do
#echo $i
#cat 'vision_rna.chr'*'.'*'.'$i'.gene_ccRE.txt' > 'vision_rna.gene_ccRE.c'$i'.txt'
#cat 'vision_rna.gene_ccRE.c'$i'.txt' | awk -F '\t' -v OFS='\t' '{if ($8==1) print $0}' > 'vision_rna.gene_ccRE.c'$i'.selected.txt'
cat 'vision_rna.gene_ccRE.c'$i'.selected.txt' | awk -F '\t' -v OFS='\t' '{if ($3=="protein_coding") print $0}' > 'vision_rna.gene_ccRE.c'$i'.selected.protein_coding.txt'
cat 'vision_rna.gene_ccRE.c'$i'.selected.txt' | awk -F '\t' -v OFS='\t' '{if ($3!="protein_coding") print $0}' > 'vision_rna.gene_ccRE.c'$i'.selected.other.txt'
wc -l 'vision_rna.gene_ccRE.c'$i'.selected.txt'
done
