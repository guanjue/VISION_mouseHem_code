#python ../snapshot_all_pre/index_set_20cell_pknorm_qda/vlookup.py -t rnaHtseqCountsall.txt -m 1 -s gencode.vM4.annotation.bed -n 4 -o rnaHtseqCountsall_withcoordinates.txt
#paste gencode.vM4.annotation.bed rnaHtseqCountsall_withcoordinates.txt > rnaHtseqCountsall_withcoordinates.txt.tmp && mv rnaHtseqCountsall_withcoordinates.txt.tmp rnaHtseqCountsall_withcoordinates.txt
#rnaHtseqCountsall.txt

cat rnaHtseqCountsall_withcoordinates.meanheader.txt > rnaHtseqCountsall_withcoordinates.mean.txt
tail -n+2 rnaHtseqCountsall_withcoordinates.txt | awk -F '\t' -v OFS='\t' '{print $1,$2,$3,$4,$5,$6,$7, ($8+$9)/2,($10+$11)/2,($12+$13)/2,($14+$15)/2,($16+$17)/2,($18+$19)/2,($20+$21)/2,($22+$23)/2,$24,($25+$26)/2,($27+$28)/2,($29+$30)/2 }' >> rnaHtseqCountsall_withcoordinates.mean.txt


Rscript get_TPM_from_HTseqCount.R rnaHtseqCountsall_withcoordinates.mean.txt rnaHtseqCountsall_withcoordinates.mean.0.txt rnaTPMall_withcoordinates.mean.0.txt

