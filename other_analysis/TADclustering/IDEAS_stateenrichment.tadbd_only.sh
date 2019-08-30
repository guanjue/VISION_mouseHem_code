### get bd bed
cat OnTADraw_pen0.1_max200_meannorm_chr*.boundary | sort -k1,1 -k2,2n > G1E.hic.boundary.bed
for i in {1..5}
do
	echo $i
	cat G1E.hic.boundary.bed | awk -F '\t' -v OFS='\t' -v level=$i '{if ($4==level && $2-10000>0) print $1, $2-10000, $3, $4, $4; else if ($4==level) print $1, 0, $3, $4, $4}' > 'G1E.hic.boundary.'$i'.bed'
done
### get tad bed
cat OnTADraw_pen0.1_max200_meannorm_chr*.bed2 | sort -k1,1 -k2,2n > G1E.hic.tad.bed
for i in {1..5}
do
	echo $i
	cat G1E.hic.tad.bed | awk -F '\t' -v OFS='\t' -v level=$i '{if ($4==level && $2>0) print $1, $2, $3, $4, $4; else if ($4==level) print $1, 0, $3, $4, $4}' > 'G1E.hic.tad.'$i'.bed'
done

### pool tad and bd
cat G1E.hic.boundary.*.bed | sort -k1,1 -k2,2n > G1E.hic.tadbdonly.bed

### loop TAD level
mkdir IDEAS_state_coverage_tadbdonly
cp G1E.hic.tadbdonly.bed G1E.hic.tadbdonly.mat.c.txt
cp G1E.hic.tadbdonly.bed G1E.hic.tadbdonly.mat.p.txt

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
bedtools coverage -a G1E.hic.tadbdonly.bed -b 'ER4.state.'$i'.merge.bed' > 'G1E.hic.tadbdonly.ER4.IDEAS.info.'$i'.txt'
### get count
cut -f7 'G1E.hic.tadbdonly.ER4.IDEAS.info.'$i'.txt' > 'G1E.hic.tadbdonly.ER4.IDEAS.region.'$i'.txt'
paste 'G1E.hic.tadbdonly.mat.c.txt' 'G1E.hic.tadbdonly.ER4.IDEAS.region.'$i'.txt' > 'G1E.hic.tadbdonly.mat.c.txt.tmp' \
&& mv 'G1E.hic.tadbdonly.mat.c.txt.tmp' 'G1E.hic.tadbdonly.mat.c.txt'
rm 'G1E.hic.tadbdonly.ER4.IDEAS.region.'$i'.txt'
### get percentage
cut -f9 'G1E.hic.tadbdonly.ER4.IDEAS.info.'$i'.txt' > 'G1E.hic.tadbdonly.ER4.IDEAS.percentage.'$i'.txt'
paste 'G1E.hic.tadbdonly.mat.p.txt' 'G1E.hic.tadbdonly.ER4.IDEAS.percentage.'$i'.txt' > 'G1E.hic.tadbdonly.mat.p.txt.tmp' \
&& mv 'G1E.hic.tadbdonly.mat.p.txt.tmp' 'G1E.hic.tadbdonly.mat.p.txt'
rm 'G1E.hic.tadbdonly.ER4.IDEAS.percentage.'$i'.txt'
#########
### mv files into folder
#mv 'ER4.state.'$i'.merge.bed' IDEAS_state_coverage_level_separate
mv 'G1E.hic.tadbdonly.ER4.IDEAS.info.'$i'.txt' IDEAS_state_coverage_tadbdonly
done


### merge states
time Rscript get_ideas_input_tadbdonly.R

### QDA clustering
time Rscript cluster.tadbdonly.qda.R

### get bb tracks
time ~/group/software/ucsc/bedToBigBed output_newid_color.bedgraph ~/group/genome/mm10/mm10.1to19_X.10kb.genome cluster.tad.qda.bb


