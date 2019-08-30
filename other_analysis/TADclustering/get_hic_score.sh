stewd('/storage/home/g/gzx103/group/HiC/mouse/processed/10kb')

cd /storage/home/g/gzx103/group/projects/3d/G1E-ER4-uninduced.merged/10kb/run_TAD_IDEAS_tad_and_bd_result

cat run_TAD_IDEAS_tad_and_bd.state | awk -F ' ' -v OFS='\t' '{if ($2=="chr19") print $2, $3, $4, $5}' > run_TAD_IDEAS_tad_and_bd.chr19.state
cat run_TAD_IDEAS_tad_and_bd.state | awk -F ' ' -v OFS='\t' '{if ($2=="chr11") print $2, $3, $4, $5}' > run_TAD_IDEAS_tad_and_bd.chr11.state
cat run_TAD_IDEAS_tad_and_bd.state | awk -F ' ' -v OFS='\t' '{print $2, $3, $4, $5}' > run_TAD_IDEAS_tad_and_bd.all.state

cat ~/group/projects/vision/mm10.10KB.bed | awk -F '\t' -v OFS='\t' '{if ($1=="chr19") print $1, $2, $3}' > mm10.10KB.chr19.bed
cat ~/group/projects/vision/mm10.10KB.bed | awk -F '\t' -v OFS='\t' '{if ($1=="chr11") print $1, $2, $3}' > mm10.10KB.chr11.bed
cat ~/group/projects/vision/mm10.10KB.bed | awk -F '\t' -v OFS='\t' '{print $1, $2, $3}' > mm10.10KB.all.bed

cat ../TAD_regions_tad_and_bd.txt | awk -F ' ' -v OFS='\t' '{print $1, $2, $3, $4}' | sort -k1,1 -k2,2n > TAD_regions_tad_and_bd.bed


R
d = read.table('mm10.10KB.chr19.bed', header=F)
d = cbind(d, c(1:dim(d)[1]))
write.table(d, 'mm10.10KB.chr19.id.bed', quote=F, col.names=F, row.names=F, sep='\t')
d = read.table('mm10.10KB.all.bed', header=F)
d = cbind(d, c(1:dim(d)[1]))
write.table(d, 'mm10.10KB.all.id.bed', quote=F, col.names=F, row.names=F, sep='\t')
d = read.table('mm10.10KB.chr11.bed', header=F)
d = cbind(d, c(1:dim(d)[1]))
write.table(d, 'mm10.10KB.chr11.id.bed', quote=F, col.names=F, row.names=F, sep='\t')


for i in {0..14}
do
echo $i
cat run_TAD_IDEAS_tad_and_bd.chr19.state | awk -F '\t' -v OFS='\t' -v state=$i '{if ($4==state) print $1, $2, $3, $4}' > 'chr19.'$i'.state.bed'
bedtools intersect -a mm10.10KB.chr19.id.bed -b 'chr19.'$i'.state.bed' -wa > 'mm10.10KB.chr19.s'$i'.bed'
done


cat /storage/home/g/gzx103/group/projects/3d/G1E-ER4-uninduced.merged/10kb/CP/G1E-ER4all3.merged.compartment.revised.sort.bedgraph \
| awk -F '\t' -v OFS='\t' -v state=$i '{if ($4<=0) print $1, $2, $3, $4}' > G1E-ER4all3.merged.CP.neg.bed
cat /storage/home/g/gzx103/group/projects/3d/G1E-ER4-uninduced.merged/10kb/CP/G1E-ER4all3.merged.compartment.revised.sort.bedgraph \
| awk -F '\t' -v OFS='\t' -v state=$i '{if ($4>0) print $1, $2, $3, $4}' > G1E-ER4all3.merged.CP.pos.bed

bedtools intersect -a mm10.10KB.chr19.id.bed -b G1E-ER4all3.merged.CP.neg.bed -wa > mm10.10KB.chr19.CPneg.bed
bedtools intersect -a mm10.10KB.chr19.id.bed -b G1E-ER4all3.merged.CP.pos.bed -wa > mm10.10KB.chr19.CPpos.bed





for i in {0..14}
do
echo $i
cat run_TAD_IDEAS_tad_and_bd.chr11.state | awk -F '\t' -v OFS='\t' -v state=$i '{if ($4==state) print $1, $2, $3, $4}' > 'chr11.'$i'.state.bed'
bedtools intersect -a mm10.10KB.chr11.id.bed -b 'chr11.'$i'.state.bed' -wa > 'mm10.10KB.chr11.s'$i'.bed'
done









for i in {0..14}
do
echo $i
cat run_TAD_IDEAS_tad_and_bd.all.state | awk -F '\t' -v OFS='\t' -v state=$i '{if ($4==state) print $1, $2, $3, $4}' > 'all.'$i'.state.bed'
bedtools intersect -a TAD_regions_tad_and_bd.bed -b 'all.'$i'.state.bed' -wa > 'TAD_regions_tad_and_bd.all.s'$i'.bed'
bedtools closest -a 'TAD_regions_tad_and_bd.all.s'$i'.bed' -b 'TAD_regions_tad_and_bd.all.s'$i'.bed' -io -d > 'TAD_regions_tad_and_bd.all.s'$i'.dist.bed'
done

