
### loop TAD level
mkdir IDEAS_state_coverage_tad_and_bd
cat ~/group/projects/vision/mm10.10KB.bed | awk -F '\t' -v OFS='\t' '{print $1, $2, $3, $1"_"$2"_"$3, $1"_"$2"_"$3}' > mm10.10kb.bed
cp mm10.10kb.bed HPC7.hic.tad_and_bd.mat.c.txt
cp mm10.10kb.bed HPC7.hic.tad_and_bd.mat.p.txt

### loop IDEAS state
for i in {0..26}
do
echo $i
### split IDEAS states
cat /storage/home/gzx103/group/projects/vision/pknorm_2_16lim_ref1mo_sample1mo_0424/get_eRP_tracks/HPC7.state.bed \
| awk -F '\t' -v OFS='\t' -v state=$i '{if ($4==state) print $0}' > 'HPC7.state.'$i'.bed'
bedtools merge -i 'HPC7.state.'$i'.bed' > 'HPC7.state.'$i'.merge.bed'
rm 'HPC7.state.'$i'.bed'
#########
bedtools coverage -a mm10.10kb.bed -b 'HPC7.state.'$i'.merge.bed' > 'HPC7.hic.tad_and_bd.HPC7.IDEAS.info.'$i'.txt'
### get count
cut -f7 'HPC7.hic.tad_and_bd.HPC7.IDEAS.info.'$i'.txt' > 'HPC7.hic.tad_and_bd.HPC7.IDEAS.region.'$i'.txt'
paste 'HPC7.hic.tad_and_bd.mat.c.txt' 'HPC7.hic.tad_and_bd.HPC7.IDEAS.region.'$i'.txt' > 'HPC7.hic.tad_and_bd.mat.c.txt.tmp' \
&& mv 'HPC7.hic.tad_and_bd.mat.c.txt.tmp' 'HPC7.hic.tad_and_bd.mat.c.txt'
rm 'HPC7.hic.tad_and_bd.HPC7.IDEAS.region.'$i'.txt'
### get percentage
cut -f9 'HPC7.hic.tad_and_bd.HPC7.IDEAS.info.'$i'.txt' > 'HPC7.hic.tad_and_bd.HPC7.IDEAS.percentage.'$i'.txt'
paste 'HPC7.hic.tad_and_bd.mat.p.txt' 'HPC7.hic.tad_and_bd.HPC7.IDEAS.percentage.'$i'.txt' > 'HPC7.hic.tad_and_bd.mat.p.txt.tmp' \
&& mv 'HPC7.hic.tad_and_bd.mat.p.txt.tmp' 'HPC7.hic.tad_and_bd.mat.p.txt'
rm 'HPC7.hic.tad_and_bd.HPC7.IDEAS.percentage.'$i'.txt'
#########
### mv files into folder
#mv 'HPC7.state.'$i'.merge.bed' IDEAS_state_coverage_level_separate
mv 'HPC7.hic.tad_and_bd.HPC7.IDEAS.info.'$i'.txt' IDEAS_state_coverage_tad_and_bd
done


cp mm10.10kb.bed ER4.hic.tad_and_bd.mat.c.txt
cp mm10.10kb.bed ER4.hic.tad_and_bd.mat.p.txt

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
bedtools coverage -a mm10.10kb.bed -b 'ER4.state.'$i'.merge.bed' > 'ER4.hic.tad_and_bd.ER4.IDEAS.info.'$i'.txt'
### get count
cut -f7 'ER4.hic.tad_and_bd.ER4.IDEAS.info.'$i'.txt' > 'ER4.hic.tad_and_bd.ER4.IDEAS.region.'$i'.txt'
paste 'ER4.hic.tad_and_bd.mat.c.txt' 'ER4.hic.tad_and_bd.ER4.IDEAS.region.'$i'.txt' > 'ER4.hic.tad_and_bd.mat.c.txt.tmp' \
&& mv 'ER4.hic.tad_and_bd.mat.c.txt.tmp' 'ER4.hic.tad_and_bd.mat.c.txt'
rm 'ER4.hic.tad_and_bd.ER4.IDEAS.region.'$i'.txt'
### get percentage
cut -f9 'ER4.hic.tad_and_bd.ER4.IDEAS.info.'$i'.txt' > 'ER4.hic.tad_and_bd.ER4.IDEAS.percentage.'$i'.txt'
paste 'ER4.hic.tad_and_bd.mat.p.txt' 'ER4.hic.tad_and_bd.ER4.IDEAS.percentage.'$i'.txt' > 'ER4.hic.tad_and_bd.mat.p.txt.tmp' \
&& mv 'ER4.hic.tad_and_bd.mat.p.txt.tmp' 'ER4.hic.tad_and_bd.mat.p.txt'
rm 'ER4.hic.tad_and_bd.ER4.IDEAS.percentage.'$i'.txt'
#########
### mv files into folder
#mv 'ER4.state.'$i'.merge.bed' IDEAS_state_coverage_level_separate
mv 'ER4.hic.tad_and_bd.ER4.IDEAS.info.'$i'.txt' IDEAS_state_coverage_tad_and_bd
done




