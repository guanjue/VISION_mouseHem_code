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
### get CP bd bed
tail -n+2 G1E-ER4all3.merged.compartment.revised.bedgraph | awk -F '\t' -v OFS='\t' '{if ($4>=0) print $1,$2,$3}' | sort -k1,1 -k2,2n > G1E.hic.CP.pos.bed
tail -n+2 G1E-ER4all3.merged.compartment.revised.bedgraph | awk -F '\t' -v OFS='\t' '{if ($4<0) print $1,$2,$3}' | sort -k1,1 -k2,2n > G1E.hic.CP.neg.bed
bedtools merge -i G1E.hic.CP.pos.bed -d 50000 > G1E.hic.CP.pos.merge.bed
bedtools merge -i G1E.hic.CP.neg.bed -d 50000 > G1E.hic.CP.neg.merge.od.bed
bedtools subtract -a G1E.hic.CP.neg.merge.od.bed -b G1E.hic.CP.pos.merge.bed > G1E.hic.CP.neg.merge.bed
cat G1E.hic.CP.pos.merge.bed G1E.hic.CP.neg.merge.bed | sort -k1,1 -k2,2n | awk -F '\t' -v OFS='\t' '{print $1,$2-10000,$2+10000}' > G1E.hic.CP.bd.bed
cat G1E.hic.CP.pos.merge.bed G1E.hic.CP.neg.merge.bed | sort -k1,1 -k2,2n | tail -1 | awk -F '\t' -v OFS='\t' '{print $1,$3-10000,$3+10000}' >> G1E.hic.CP.bd.bed
### get CP bed
cat G1E.hic.CP.pos.merge.bed G1E.hic.CP.neg.merge.bed | sort -k1,1 -k2,2n > G1E.hic.CP.pos_neg.merge.bed
bedtools subtract -a G1E.hic.CP.pos_neg.merge.bed -b G1E.hic.CP.bd.bed > G1E.hic.CP.pos_neg.merge.nobd.bed


### get CP & CP bd
cat G1E.hic.CP.pos_neg.merge.nobd.bed G1E.hic.CP.bd.bed | sort -u | sort -k1,1 -k2,2n > G1E.hic.CP.CPbd.bed

### split tads by tad levels
for i in {1..5}
do
	echo $i
	cat G1E.hic.tad.bed | awk -F '\t' -v OFS='\t' -v level=$i '{if ($4==level && $2>0) print $1, $2, $3, $4, $4; else if ($4==level) print $1, 0, $3, $4, $4}' > 'G1E.hic.tad.'$i'.bed'
done

### pool tad and bd
rm G1E.hic.tad.*.sub.bed
cat G1E.hic.tad.*.bed G1E.hic.boundary.*.bed | sort -u | sort -k1,1 -k2,2n > G1E.hic.tad_and_bd.od.bed

### get CP CP bd uniques
bedtools subtract -a G1E.hic.CP.CPbd.bed -b G1E.hic.tad_and_bd.od.bed > G1E.hic.CP.CPbd.OK.bed
bedtools subtract -a G1E.hic.CP.pos_neg.merge.bed -b G1E.hic.tad_and_bd.od.bed > G1E.hic.CP.notad.OK.bed


### remove low level tads regions
for i in {1..5}
do
	cp 'G1E.hic.tad.'$i'.bed' 'G1E.hic.tad.'$i'.sub.bed'
	for j in {1..5}
	do
		if (( i < j ))
		then
			echo $i $j
			bedtools subtract -a 'G1E.hic.tad.'$i'.sub.bed' -b 'G1E.hic.tad.'$j'.bed' > 'G1E.hic.tad.'$i'.sub.bed.tmp'
			mv 'G1E.hic.tad.'$i'.sub.bed.tmp' 'G1E.hic.tad.'$i'.sub.bed'
		fi
	done
	### loop boundary
	for k in {1..5}
	do
		echo $i $k
		bedtools subtract -a 'G1E.hic.tad.'$i'.sub.bed' -b 'G1E.hic.boundary.'$k'.bed' > 'G1E.hic.tad.'$i'.sub.bed.tmp'
		mv 'G1E.hic.tad.'$i'.sub.bed.tmp' 'G1E.hic.tad.'$i'.sub.bed'
	done
done

### pool tad and bd
cat G1E.hic.tad.*.sub.bed G1E.hic.boundary.*.bed | sort -k1,1 -k2,2n > G1E.hic.tad_and_bd.noCP.bed

### get CP split tad bd
bedtools intersect -a G1E.hic.tad_and_bd.noCP.bed -b G1E.hic.CP.CPbd.bed -c > G1E.hic.tad_and_bd.CPcount.bed
### get no CP split tad & tad bd
cat G1E.hic.tad_and_bd.CPcount.bed | awk -F '\t' -v OFS='\t' -v level=$i '{if ($6==1) print $1, $2, $3, $4, $5}' > G1E.hic.tad_and_bd.OK.bed
### split tad and tad bd by CP
cat G1E.hic.tad_and_bd.CPcount.bed | awk -F '\t' -v OFS='\t' -v level=$i '{if ($6>1) print $1, $2, $3, $4, $5}' > G1E.hic.tad_and_bd.CPsplit.0.bed
bedtools intersect -a G1E.hic.CP.CPbd.bed -b G1E.hic.tad_and_bd.CPsplit.bed > G1E.hic.tad_and_bd.CPsplit.1.bed
bedtools subtract -a G1E.hic.tad_and_bd.CPsplit.bed -b G1E.hic.tad_and_bd.CPsplit.1.bed > G1E.hic.tad_and_bd.CPsplit.2.bed
cat G1E.hic.tad_and_bd.CPsplit.1.bed G1E.hic.tad_and_bd.CPsplit.2.bed | sort -u | sort -k1,1 -k2,2n > G1E.hic.tad_and_bd.CPsplit.OK.bed

### get all split region bed file
#cat G1E.hic.tad_and_bd.OK.bed G1E.hic.CP.CPbd.OK.bed G1E.hic.tad_and_bd.CPsplit.OK.bed | sort -u | sort -k1,1 -k2,2n \
#| awk -F '\t' -v OFS='\t' '{print $1, $2, $3, 1, 1}' > G1E.hic.tad_and_bd_and_CP.bed
cat G1E.hic.tad_and_bd.noCP.bed G1E.hic.CP.notad.OK.bed | sort -u | sort -k1,1 -k2,2n \
| awk -F '\t' -v OFS='\t' '{print $1, $2, $3, 1, 1}' > G1E.hic.tad_and_bd_and_CP.bed



### loop TAD level
mkdir IDEAS_state_coverage_tad_and_bd_and_CP
cp G1E.hic.tad_and_bd_and_CP.bed G1E.hic.tad_and_bd_and_CP.mat.c.txt
cp G1E.hic.tad_and_bd_and_CP.bed G1E.hic.tad_and_bd_and_CP.mat.p.txt

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
bedtools coverage -a G1E.hic.tad_and_bd_and_CP.bed -b 'ER4.state.'$i'.merge.bed' > 'G1E.hic.tad_and_bd_and_CP.ER4.IDEAS.info.'$i'.txt'
### get count
cut -f7 'G1E.hic.tad_and_bd_and_CP.ER4.IDEAS.info.'$i'.txt' > 'G1E.hic.tad_and_bd_and_CP.ER4.IDEAS.region.'$i'.txt'
paste 'G1E.hic.tad_and_bd_and_CP.mat.c.txt' 'G1E.hic.tad_and_bd_and_CP.ER4.IDEAS.region.'$i'.txt' > 'G1E.hic.tad_and_bd_and_CP.mat.c.txt.tmp' \
&& mv 'G1E.hic.tad_and_bd_and_CP.mat.c.txt.tmp' 'G1E.hic.tad_and_bd_and_CP.mat.c.txt'
rm 'G1E.hic.tad_and_bd_and_CP.ER4.IDEAS.region.'$i'.txt'
### get percentage
cut -f9 'G1E.hic.tad_and_bd_and_CP.ER4.IDEAS.info.'$i'.txt' > 'G1E.hic.tad_and_bd_and_CP.ER4.IDEAS.percentage.'$i'.txt'
paste 'G1E.hic.tad_and_bd_and_CP.mat.p.txt' 'G1E.hic.tad_and_bd_and_CP.ER4.IDEAS.percentage.'$i'.txt' > 'G1E.hic.tad_and_bd_and_CP.mat.p.txt.tmp' \
&& mv 'G1E.hic.tad_and_bd_and_CP.mat.p.txt.tmp' 'G1E.hic.tad_and_bd_and_CP.mat.p.txt'
rm 'G1E.hic.tad_and_bd_and_CP.ER4.IDEAS.percentage.'$i'.txt'
#########
### mv files into folder
#mv 'ER4.state.'$i'.merge.bed' IDEAS_state_coverage_level_separate
mv 'G1E.hic.tad_and_bd_and_CP.ER4.IDEAS.info.'$i'.txt' IDEAS_state_coverage_tad_and_bd_and_CP
done



255,255,255
148,148,148
18,13,197
40,40,215
204.2,33,151.4
22,138,22
170,179.5,0
126.5,175.5,65
203.454545454545,100.636363636364,49.5454545454545

255,255,255
148,148,148
18,13,197
157.285714285714,35,169.571428571429
22,138,22
170,179.5,0
191.615384615385,112.153846153846,51.9230769230769

255,255,255
148,148,148
18,13,197
157.285714285714,35,169.571428571429
22,138,22
188.733333333333,121.133333333333,45



255,255,255
86.8333333333333,149.166666666667,48.3333333333333
148,148,148
29,26.5,206
232.25,166,73.75
197.333333333333,40,218.666666666667
229,55.125,26.625


255,255,255
148,148,148
18,13,197
153,40.8333333333333,184.333333333333
22,138,22
170,179.5,0
191,104.142857142857,54


