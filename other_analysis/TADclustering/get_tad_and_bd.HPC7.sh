### get bd bed
cat OnTADraw_pen0.1_max200_meannorm_chr*.boundary | sort -k1,1 -k2,2n > HPC7.hic.boundary.bed
for i in {1..5}
do
	echo $i
	cat HPC7.hic.boundary.bed | awk -F '\t' -v OFS='\t' -v level=$i '{if ($4==level && $2-10000>0) print $1, $2-10000, $3, $4, $4; else if ($4==level) print $1, 0, $3, $4, $4}' > 'HPC7.hic.boundary.'$i'.bed'
done
### get tad bed
cat OnTADraw_pen0.1_max200_meannorm_chr*.bed2 | sort -k1,1 -k2,2n > HPC7.hic.tad.bed
for i in {1..5}
do
	echo $i
	cat HPC7.hic.tad.bed | awk -F '\t' -v OFS='\t' -v level=$i '{if ($4==level && $2>0) print $1, $2, $3, $4, $4; else if ($4==level) print $1, 0, $3, $4, $4}' > 'HPC7.hic.tad.'$i'.bed'
done

### pool tad and bd
rm HPC7.hic.tad.*.sub.bed
cat HPC7.hic.*.*.bed | sort -k1,1 -k2,2n > HPC7.hic.tad_and_bd.bed

### remove low level tads regions
for i in {1..5}
do
cp 'HPC7.hic.tad.'$i'.bed' 'HPC7.hic.tad.'$i'.sub.bed'
for j in {1..5}
do	
if (( i < j ))
then
echo $i $j
bedtools subtract -a 'HPC7.hic.tad.'$i'.sub.bed' -b 'HPC7.hic.tad.'$j'.bed' > 'HPC7.hic.tad.'$i'.sub.bed.tmp'
mv 'HPC7.hic.tad.'$i'.sub.bed.tmp' 'HPC7.hic.tad.'$i'.sub.bed'
fi
done
### loop boundary
for k in {1..5}
do	
echo $i $k
bedtools subtract -a 'HPC7.hic.tad.'$i'.sub.bed' -b 'HPC7.hic.boundary.'$k'.bed' > 'HPC7.hic.tad.'$i'.sub.bed.tmp'
mv 'HPC7.hic.tad.'$i'.sub.bed.tmp' 'HPC7.hic.tad.'$i'.sub.bed'
done
done

### pool tad and bd
cat HPC7.hic.tad.*.sub.bed HPC7.hic.boundary.*.bed | sort -k1,1 -k2,2n > HPC7.hic.tad_and_bd.bed

### loop TAD level
mkdir IDEAS_state_coverage_tad_and_bd
cp HPC7.hic.tad_and_bd.bed HPC7.hic.tad_and_bd.mat.c.txt
cp HPC7.hic.tad_and_bd.bed HPC7.hic.tad_and_bd.mat.p.txt

### loop IDEAS state
for i in {0..26}
do
echo $i
### split IDEAS states
#cat /storage/home/gzx103/group/projects/vision/pknorm_2_16lim_ref1mo_sample1mo_0424/get_eRP_tracks/ER4.state.bed \
#| awk -F '\t' -v OFS='\t' -v state=$i '{if ($4==state) print $0}' > 'ER4.state.'$i'.bed'
#bedtools merge -i 'ER4.state.'$i'.bed' > 'ER4.state.'$i'.merge.bed'
#rm 'ER4.state.'$i'.bed'
#########
bedtools coverage -a HPC7.hic.tad_and_bd.bed -b 'ER4.state.'$i'.merge.bed' > 'HPC7.hic.tad_and_bd.ER4.IDEAS.info.'$i'.txt'
### get count
cut -f7 'HPC7.hic.tad_and_bd.ER4.IDEAS.info.'$i'.txt' > 'HPC7.hic.tad_and_bd.ER4.IDEAS.region.'$i'.txt'
paste 'HPC7.hic.tad_and_bd.mat.c.txt' 'HPC7.hic.tad_and_bd.ER4.IDEAS.region.'$i'.txt' > 'HPC7.hic.tad_and_bd.mat.c.txt.tmp' \
&& mv 'HPC7.hic.tad_and_bd.mat.c.txt.tmp' 'HPC7.hic.tad_and_bd.mat.c.txt'
rm 'HPC7.hic.tad_and_bd.ER4.IDEAS.region.'$i'.txt'
### get percentage
cut -f9 'HPC7.hic.tad_and_bd.ER4.IDEAS.info.'$i'.txt' > 'HPC7.hic.tad_and_bd.ER4.IDEAS.percentage.'$i'.txt'
paste 'HPC7.hic.tad_and_bd.mat.p.txt' 'HPC7.hic.tad_and_bd.ER4.IDEAS.percentage.'$i'.txt' > 'HPC7.hic.tad_and_bd.mat.p.txt.tmp' \
&& mv 'HPC7.hic.tad_and_bd.mat.p.txt.tmp' 'HPC7.hic.tad_and_bd.mat.p.txt'
rm 'HPC7.hic.tad_and_bd.ER4.IDEAS.percentage.'$i'.txt'
#########
### mv files into folder
#mv 'ER4.state.'$i'.merge.bed' IDEAS_state_coverage_level_separate
mv 'HPC7.hic.tad_and_bd.ER4.IDEAS.info.'$i'.txt' IDEAS_state_coverage_tad_and_bd
done


