#PBS -l nodes=1:ppn=8
#PBS -l walltime=20:00:00
#PBS -j oe
#PBS -A yzz2_e_g_sc_default
#PBS -l pmem=16gb

module load gcc/5.3.1
module load python/2.7.14-anaconda5.0.1
module load bedtools

cd /storage/home/g/gzx103/group/projects/vision/rna/gene_ccRE_200bp_to_ccRE
rm all_converted.sort.g1234.uniq.txt
cat all_converted.sort.g*.uniq.txt > all_converted.sort.g1234.uniq.txt

mkdir gene_ccREs_bw_1216
while read LINE
do
	outputname=$(echo $LINE | awk -F ':' '{print $1"."$2}')
	echo $outputname
	###
	#cat vision_rna.gene_ccRE.selected.all.reprod_count120718_withNames.txt | awk -F '\t' -v OFS='\t' -v g=$LINE '{if ($3==g) print $0}' | awk -F '\t' -v OFS='\t' '{ split($10, chars, ""); num=0; for (i=1; i <= length($10); i++) {num = num + int(chars[i])}; $6=sprintf("%d",$6); print $1, $6-200, $6, $8, num}' | sort -k1,1 -k2,2n > $outputname'.cor_repr.txt'
	cat all_converted.sort.g1234.uniq.txt | awk -F '\t' -v OFS='\t' -v g=$LINE '{if ($3==g) print $1, $4, $5, $7, $8}' | sort -k1,1 -k2,2n  > $outputname'.ccRE.allinfo.txt'
	time python extract_bedgraph_ccRE.py -i $outputname'.ccRE.allinfo.txt' -o $outputname'.ccRE'
	~/group/software/ucsc/bedGraphToBigWig $outputname'.ccRE.cor.bedgraph' ~/group/projects/vision/input_norm/mm10.chrom.sizes $outputname'.ccRE.cor.bw'
	~/group/software/ucsc/bedGraphToBigWig $outputname'.ccRE.repr.bedgraph' ~/group/projects/vision/input_norm/mm10.chrom.sizes $outputname'.ccRE.repr.bw'
	mv $outputname'.ccRE.repr.bw' gene_ccREs_bw_1216
	mv $outputname'.ccRE.repr.bedgraph' gene_ccREs_bw_1216
	mv $outputname'.ccRE.allinfo.txt' gene_ccREs_bw_1216
	mv $outputname'.ccRE.cor.bw' gene_ccREs_bw_1216
	mv $outputname'.ccRE.cor.bedgraph' gene_ccREs_bw_1216
#done < vision_rna.gene_ccRE.uniqNames.test.txt
done < vision_rna.gene_ccRE.uniqNames.1207.txt
