time Rscript get_tad_score.R

### split TADs
for i in {1..5}
do
	echo $i
	cat bed_tads_bed_withscore.bed | awk -F '\t' -v OFS='\t' -v level=$i '{if ($11==level) print $0}' > 'bed_tads_bed_withscore.'$i'.bed'
done

### get IDEAS bed
mkdir IDEAS_state_coverage_level_separate
sort -k1,1 -k2,2n bed_tads_bed_withscore.bed | awk -F '\t' -v OFS='\t' '{print $1, $2, $3, $10, $11}' > bed_tads_bed_withscore.sort.bed 
cp bed_tads_bed_withscore.sort.bed tad.IDEAS.ER4.mat.c.txt
cp bed_tads_bed_withscore.sort.bed tad.IDEAS.ER4.mat.p.txt

### get different TAD level peaks
for j in {1..5}
do
echo $j
cat bed_tads_bed_withscore.sort.bed | awk -F '\t' -v OFS='\t' -v level=$j '{if ($5==level) print $0}' > 'tad.level.'$j'.bed'
done

### remove low level tads regions
for i in {1..5}
do
cp 'tad.level.'$i'.bed' 'tad.level.'$i'.sub.bed'
for j in {1..5}
do	
if (( i < j ))
then
echo $i $j
bedtools subtract -a 'tad.level.'$i'.sub.bed' -b 'tad.level.'$j'.bed' > 'tad.level.'$i'.sub.bed.tmp'
mv 'tad.level.'$i'.sub.bed.tmp' 'tad.level.'$i'.sub.bed'
fi
done
done

### loop TAD level
for j in {1..5}
do
echo $j
cp 'tad.level.'$j'.sub.bed' 'tad.level.'$j'.mat.c.txt'
cp 'tad.level.'$j'.sub.bed' 'tad.level.'$j'.mat.p.txt'
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
bedtools coverage -a 'tad.level.'$j'.sub.bed' -b 'ER4.state.'$i'.merge.bed' > 'tad.'$j'.ER4.IDEAS.info.'$i'.txt'
### get count
cut -f7 'tad.'$j'.ER4.IDEAS.info.'$i'.txt' > 'tad.'$j'.ER4.IDEAS.region.'$i'.txt'
paste 'tad.level.'$j'.mat.c.txt' 'tad.'$j'.ER4.IDEAS.region.'$i'.txt' > 'tad.level.'$j'.mat.c.txt.tmp' \
&& mv 'tad.level.'$j'.mat.c.txt.tmp' 'tad.level.'$j'.mat.c.txt'
rm 'tad.'$j'.ER4.IDEAS.region.'$i'.txt'
### get percentage
cut -f9 'tad.'$j'.ER4.IDEAS.info.'$i'.txt' > 'tad.'$j'.ER4.IDEAS.percentage.'$i'.txt'
paste 'tad.level.'$j'.mat.p.txt' 'tad.'$j'.ER4.IDEAS.percentage.'$i'.txt' > 'tad.level.'$j'.mat.p.txt.tmp' \
&& mv 'tad.level.'$j'.mat.p.txt.tmp' 'tad.level.'$j'.mat.p.txt'
rm 'tad.'$j'.ER4.IDEAS.percentage.'$i'.txt'
#########
### mv files into folder
#mv 'ER4.state.'$i'.merge.bed' IDEAS_state_coverage_level_separate
mv 'tad.'$j'.ER4.IDEAS.info.'$i'.txt' IDEAS_state_coverage_level_separate
done
done



### loop TAD level
for j in {1..5}
do
echo $j
cp 'tad.level.'$j'.bed' 'tad.level.'$j'.full.mat.c.txt'
cp 'tad.level.'$j'.bed' 'tad.level.'$j'.full.mat.p.txt'
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
bedtools coverage -a 'tad.level.'$j'.bed' -b 'ER4.state.'$i'.merge.bed' > 'tad.'$j'.ER4.IDEAS.info.'$i'.full.txt'
### get count
cut -f7 'tad.'$j'.ER4.IDEAS.info.'$i'.full.txt' > 'tad.'$j'.ER4.IDEAS.region.'$i'.full.txt'
paste 'tad.level.'$j'.full.mat.c.txt' 'tad.'$j'.ER4.IDEAS.region.'$i'.full.txt' > 'tad.level.'$j'.full.mat.c.txt.tmp' \
&& mv 'tad.level.'$j'.full.mat.c.txt.tmp' 'tad.level.'$j'.full.mat.c.txt'
rm 'tad.'$j'.ER4.IDEAS.region.'$i'.full.txt'
### get percentage
cut -f9 'tad.'$j'.ER4.IDEAS.info.'$i'.full.txt' > 'tad.'$j'.ER4.IDEAS.percentage.'$i'.full.txt'
paste 'tad.level.'$j'.full.mat.p.txt' 'tad.'$j'.ER4.IDEAS.percentage.'$i'.full.txt' > 'tad.level.'$j'.full.mat.p.txt.tmp' \
&& mv 'tad.level.'$j'.full.mat.p.txt.tmp' 'tad.level.'$j'.full.mat.p.txt'
rm 'tad.'$j'.ER4.IDEAS.percentage.'$i'.full.txt'
#########
### mv files into folder
#mv 'ER4.state.'$i'.merge.bed' IDEAS_state_coverage_level_separate
mv 'tad.'$j'.ER4.IDEAS.info.'$i'.full.txt' IDEAS_state_coverage_level_separate
done
done














### get IDEAS bed
mkdir IDEAS_state_coverage
sort -k1,1 -k2,2n bed_tads_bed_withscore.bed | awk -F '\t' -v OFS='\t' '{print $1, $2, $3, $10, $11}' > bed_tads_bed_withscore.sort.bed 
cp bed_tads_bed_withscore.sort.bed tad.IDEAS.ER4.mat.c.txt
cp bed_tads_bed_withscore.sort.bed tad.IDEAS.ER4.mat.p.txt
### loop IDEAS state
for i in {0..26}
do
echo $i
### split IDEAS states
cat /storage/home/gzx103/group/projects/vision/pknorm_2_16lim_ref1mo_sample1mo_0424/get_eRP_tracks/ER4.state.bed \
| awk -F '\t' -v OFS='\t' -v state=$i '{if ($4==state) print $0}' > 'ER4.state.'$i'.bed'
bedtools merge -i 'ER4.state.'$i'.bed' > 'ER4.state.'$i'.merge.bed'
rm 'ER4.state.'$i'.bed'
#########
bedtools coverage -a bed_tads_bed_withscore.sort.bed -b 'ER4.state.'$i'.merge.bed' > 'tad_coverage.ER4.IDEAS.info.'$i'.txt'
cut -f7 'tad_coverage.ER4.IDEAS.info.'$i'.txt' > 'tad_coverage.ER4.IDEAS.region.'$i'.txt'
cut -f9 'tad_coverage.ER4.IDEAS.info.'$i'.txt' > 'tad_coverage.ER4.IDEAS.percentage.'$i'.txt'
paste tad.IDEAS.ER4.mat.c.txt 'tad_coverage.ER4.IDEAS.region.'$i'.txt' > tad.IDEAS.ER4.mat.c.txt.tmp \
&& mv tad.IDEAS.ER4.mat.c.txt.tmp tad.IDEAS.ER4.mat.c.txt
paste tad.IDEAS.ER4.mat.p.txt 'tad_coverage.ER4.IDEAS.percentage.'$i'.txt' > tad.IDEAS.ER4.mat.p.txt.tmp \
&& mv tad.IDEAS.ER4.mat.p.txt.tmp tad.IDEAS.ER4.mat.p.txt
rm 'tad_coverage.ER4.IDEAS.region.'$i'.txt'
rm 'tad_coverage.ER4.IDEAS.percentage.'$i'.txt'
#########
### mv files into folder
#cp 'ER4.state.'$i'.merge.bed' IDEAS_state_coverage
mv 'tad_coverage.ER4.IDEAS.info.'$i'.txt' IDEAS_state_coverage
done






