### get reads count
while read LINE
do
	name=$(echo "$LINE" | awk -F '\t' -v OFS='\t' '{print $1}'  | awk -F '.' -v OFS='.' '{print $1,$2}')
	echo $name
	samtools sort $name.bam $name.sorted
	echo 'done with sorting'
	samtools index $name.sorted.bam
	echo 'done with indexing'
	bedtools bamtobed -i $name.sorted.bam > $name.bambed
	echo 'done with bamtobed'
	cat $name.bambed | awk -F '\t' -v OFS='\t' '{if ($6=="+") print $1,$2,$2+1,$4,$5,$6; else print $1,$3-1,$3,$4,$5,$6}' > $name.bam.5p.bed
	sort -k1,1 -k2,2n $name.bam.5p.bed > $name.bam.5p.sort.bed
	bedtools intersect -a /storage/home/lua137/work/VISION/200_noblack.bed -b $name.bam.5p.sort.bed -wa -c > $name.bamtobed5endintersect
done < info_table_input.txt


### total mean normalization
while read LINE
do
	sig1=$(echo "$LINE" | awk '{print $1}')
	sig2=$(echo "$LINE" | awk '{print $2}')
	echo $sig1
	echo $sig2
	time Rscript normalize_input.R $sig1 $sig2 $sig1'.totalmean.txt'
done < info_table_input.txt

### calculate the mean
paste *'.totalmean.txt' | awk -F '\t' -v OFS='\t' '{sum=0; for(i=1; i<=NF; i++){sum+=$i}; sum/=NF; print sum}' > merged_normed_input_mean.txt

### round the mean
cat merged_normed_input_mean.txt | awk -F '\t' -v OFS='\t' '{printf "%1.0f\n", $1}' > merged_normed_input_mean.rounding.txt

