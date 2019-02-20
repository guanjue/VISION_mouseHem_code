sort -k1,1 -k2,2n rnaTPM.2KB.txt > rnaTPM.2KB.sort.txt
bedtools intersect -a rnaTPM.2KB.sort.txt -b vision_cres.bed -wao > rnaTPM.2KB.ccRE.txt


sort -k1,1 -k2,2n rnaTPM.2MB.txt > rnaTPM.2MB.sort.txt
bedtools intersect -a rnaTPM.2MB.sort.txt -b vision_cres.bed -wao > rnaTPM.2MB.ccRE.txt





