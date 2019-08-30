### get 10kb bins
tail -n+2 G1E-ER4all3.merged.compartment.revised.bedgraph | awk -F '\t' -v OFS='\t' '{print $1, $2, $3, $1"_"$2"_"$3}' | sort -k1,1 -k2,2n > G1E_ER4.CP.bed

### get CP bw
tail -n+2 G1E-ER4all3.merged.compartment.revised.bedgraph | sort -k1,1 -k2,2n > G1E-ER4all3.merged.compartment.revised.sort.bedgraph
~/group/software/ucsc/bedGraphToBigWig G1E-ER4all3.merged.compartment.revised.sort.bedgraph ~/group/genome/mm10/mm10.1to19_X.genome G1E-ER4all3.CP.bw

### get 10kb bins CP scores
~/group/software/ucsc/bigWigAverageOverBed G1E-ER4all3.CP.bw G1E_ER4.CP.bed G1E_ER4.CP.bed.CP.tab
cat G1E_ER4.CP.bed.CP.tab | cut -f1 | awk -F '_' -v OFS='\t' '{print $1, $2, $3}' > G1E_ER4.CP.bed.CP.tab.bed
paste G1E_ER4.CP.bed.CP.tab.bed G1E_ER4.CP.bed.CP.tab | sort -k1,1 -k2,2n | cut -f9 | awk -F '\t' -v OFS='\t' '{if ($1>0) print ($1-0.6068174)/0.2895379*1.855611+1.019424; else print 0}'> G1E_ER4.CP.bed.CP.pos.tab.txt
paste G1E_ER4.CP.bed.CP.tab.bed G1E_ER4.CP.bed.CP.tab | sort -k1,1 -k2,2n | cut -f9 | awk -F '\t' -v OFS='\t' '{if ($1<0) print ((-$1)-0.6688068)/0.2605836*1.855611+1.019424; else print 0}'> G1E_ER4.CP.bed.CP.neg.tab.txt

### Get IDEAS states counts
mkdir IDEAS_state_coverage_10kb
cat G1E_ER4.CP.bed | awk -F '\t' -v OFS='\t' '{print $1, $2, $3, $4, $4}' > G1E_ER4.10kb.mat.c.txt
cat G1E_ER4.CP.bed | awk -F '\t' -v OFS='\t' '{print $1, $2, $3, $4, $4}' > G1E_ER4.10kb.mat.p.txt
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
bedtools coverage -a G1E_ER4.CP.bed -b '../ER4.state.'$i'.merge.bed' > 'G1E_ER4.IDEAS.info.'$i'.txt'
### get count
cut -f6 'G1E_ER4.IDEAS.info.'$i'.txt' > 'G1E_ER4.IDEAS.region.'$i'.txt'
paste G1E_ER4.10kb.mat.c.txt 'G1E_ER4.IDEAS.region.'$i'.txt' > G1E_ER4.10kb.mat.c.txt.tmp \
&& mv G1E_ER4.10kb.mat.c.txt.tmp G1E_ER4.10kb.mat.c.txt
rm 'G1E_ER4.IDEAS.region.'$i'.txt'
### get percentage
cut -f8 'G1E_ER4.IDEAS.info.'$i'.txt' > 'G1E_ER4.IDEAS.percentage.'$i'.txt'
paste G1E_ER4.10kb.mat.p.txt 'G1E_ER4.IDEAS.percentage.'$i'.txt' > G1E_ER4.10kb.mat.p.txt.tmp \
&& mv G1E_ER4.10kb.mat.p.txt.tmp G1E_ER4.10kb.mat.p.txt
rm 'G1E_ER4.IDEAS.percentage.'$i'.txt'
#########
### mv files into folder
#mv 'ER4.state.'$i'.merge.bed' IDEAS_state_coverage_level_separate
mv 'G1E_ER4.IDEAS.info.'$i'.txt' IDEAS_state_coverage_10kb
done

