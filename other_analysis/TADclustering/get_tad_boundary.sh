cat OnTADraw_pen0.1_max200_meannorm_chr*.boundary | sort -k1,1 -k2,2n > G1E.tad.boundary.bed

for i in {1..5}
do
	echo $i
	cat G1E.tad.boundary.bed | awk -F '\t' -v OFS='\t' -v level=$i '{if ($4==level && $2>0) print $1, $2, $3, $4, $4; else if ($4==level) print $1, 0, $3, $4, $4}' > 'G1E.tad.boundary.'$i'.bed'
done


### loop TAD level
mkdir IDEAS_state_coverage_bd
for j in {1..5}
do
echo $j
cp 'G1E.tad.boundary.'$j'.bed' 'tad.bd.'$j'.mat.c.txt'
cp 'G1E.tad.boundary.'$j'.bed' 'tad.bd.'$j'.mat.p.txt'
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
bedtools coverage -a 'G1E.tad.boundary.'$j'.bed' -b 'ER4.state.'$i'.merge.bed' > 'tad.bd.'$j'.ER4.IDEAS.info.'$i'.txt'
### get count
cut -f7 'tad.bd.'$j'.ER4.IDEAS.info.'$i'.txt' > 'tad.bd.'$j'.ER4.IDEAS.region.'$i'.txt'
paste 'tad.bd.'$j'.mat.c.txt' 'tad.bd.'$j'.ER4.IDEAS.region.'$i'.txt' > 'tad.bd.'$j'.mat.c.txt.tmp' \
&& mv 'tad.bd.'$j'.mat.c.txt.tmp' 'tad.bd.'$j'.mat.c.txt'
rm 'tad.bd.'$j'.ER4.IDEAS.region.'$i'.txt'
### get percentage
cut -f9 'tad.bd.'$j'.ER4.IDEAS.info.'$i'.txt' > 'tad.bd.'$j'.ER4.IDEAS.percentage.'$i'.txt'
paste 'tad.bd.'$j'.mat.p.txt' 'tad.bd.'$j'.ER4.IDEAS.percentage.'$i'.txt' > 'tad.bd.'$j'.mat.p.txt.tmp' \
&& mv 'tad.bd.'$j'.mat.p.txt.tmp' 'tad.bd.'$j'.mat.p.txt'
rm 'tad.bd.'$j'.ER4.IDEAS.percentage.'$i'.txt'
#########
### mv files into folder
#mv 'ER4.state.'$i'.merge.bed' IDEAS_state_coverage_level_separate
mv 'tad.bd.'$j'.ER4.IDEAS.info.'$i'.txt' IDEAS_state_coverage_bd
done
done


time Rscript get_ideas_input_bd.R

time bash run_IDEAS_bd.sh

time bash get_TAD_graph.sh


