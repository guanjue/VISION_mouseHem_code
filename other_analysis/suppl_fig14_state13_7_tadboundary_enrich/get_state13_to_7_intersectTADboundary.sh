### get 200bp bins
tail -n+2 ../pknorm_2_16lim_ref1mo_sample1mo_0424/pknorm_2_16lim_ref1mo_0424_lesshet.state | awk -F ' ' -v OFS='\t' '{print $2,$3,$4,$10}' > ER4_state.200bpbin.bed
tail -n+2 ../pknorm_2_16lim_ref1mo_sample1mo_0424/pknorm_2_16lim_ref1mo_0424_lesshet.state | awk -F ' ' -v OFS='\t' '{print $2,$3,$4,$15}' > HPC7_state.200bpbin.bed
tail -n+2 ../pknorm_2_16lim_ref1mo_sample1mo_0424/pknorm_2_16lim_ref1mo_0424_lesshet.state | awk -F ' ' -v OFS='\t' '{print $2,$3,$4,$15,$10}' > HPC7_ER4_state.200bpbin.bed

### state 13 to 7 pks
cat HPC7_ER4_state.200bpbin.bed | awk -F '\t' -v OFS='\t' '{if ($4==13 && $5==13) print $1,$2,$3,$4"_"$5}' > HPC7.13.ER4.13.200bpbin.bed
cat HPC7_ER4_state.200bpbin.bed | awk -F '\t' -v OFS='\t' '{if ($4==13 && $5==7) print $1,$2,$3,$4"_"$5}' > HPC7.13.ER4.7.200bpbin.bed
bedtools merge -i HPC7.13.ER4.13.200bpbin.bed > HPC7.13.ER4.13.mergedpk.bed
bedtools merge -i HPC7.13.ER4.7.200bpbin.bed > HPC7.13.ER4.7.mergedpk.bed

### TAD boundary
cat HPC7.OnTad.joint.boundary.bed ER4.OnTad.joint.boundary.bed | awk -F '\t' -v OFS='\t' '{print $1,$2,$3}' | sort -u | sort -k1,1 -k2,2n > HPC7.ER4.uniq.OnTad.joint.boundary.bed
bedtools intersect -a HPC7.ER4.uniq.OnTad.joint.boundary.bed -b HPC7.OnTad.joint.boundary.bed -c > HPC7.ER4.uniq.OnTad.joint.boundary.HPC7count.bed
bedtools intersect -a HPC7.ER4.uniq.OnTad.joint.boundary.bed -b ER4.OnTad.joint.boundary.bed -c > HPC7.ER4.uniq.OnTad.joint.boundary.ER4count.bed
paste HPC7.ER4.uniq.OnTad.joint.boundary.HPC7count.bed HPC7.ER4.uniq.OnTad.joint.boundary.ER4count.bed | awk -F '\t' -v OFS='\t' '{if ($4!=0 && $8==0) print $1,$2,$3}' > HPC7_1.ER4_0.uniq.OnTad.joint.boundary.bed
paste HPC7.ER4.uniq.OnTad.joint.boundary.HPC7count.bed HPC7.ER4.uniq.OnTad.joint.boundary.ER4count.bed | awk -F '\t' -v OFS='\t' '{if ($4!=0 && $8!=0) print $1,$2,$3}' > HPC7_1.ER4_1.uniq.OnTad.joint.boundary.bed
paste HPC7.ER4.uniq.OnTad.joint.boundary.HPC7count.bed HPC7.ER4.uniq.OnTad.joint.boundary.ER4count.bed | awk -F '\t' -v OFS='\t' '{if ($4==0 && $8!=0) print $1,$2,$3}' > HPC7_0.ER4_1.uniq.OnTad.joint.boundary.bed

### dif bins

bedtools intersect -a HPC7_1.ER4_0.uniq.OnTad.joint.boundary.bed -b HPC7.13.ER4.13.mergedpk.bed -c > HPC7_1.ER4_0.HPC7_13.ER4_13.uniq.OnTad.joint.boundary.bed
bedtools intersect -a HPC7_1.ER4_0.uniq.OnTad.joint.boundary.bed -b HPC7.13.ER4.7.mergedpk.bed -c > HPC7_1.ER4_0.HPC7_13.ER4_7.uniq.OnTad.joint.boundary.bed
bedtools intersect -a HPC7_1.ER4_1.uniq.OnTad.joint.boundary.bed -b HPC7.13.ER4.13.mergedpk.bed -c > HPC7_1.ER4_1.HPC7_13.ER4_13.uniq.OnTad.joint.boundary.bed
bedtools intersect -a HPC7_1.ER4_1.uniq.OnTad.joint.boundary.bed -b HPC7.13.ER4.7.mergedpk.bed -c > HPC7_1.ER4_1.HPC7_13.ER4_7.uniq.OnTad.joint.boundary.bed

exp_win=(0 10000 20000 30000 40000 50000 60000 70000 80000 90000 100000)

for i in ${exp_win[@]}; 
do
echo $i
cat HPC7_1.ER4_0.uniq.OnTad.joint.boundary.bed | awk -F '\t' -v OFS='\t' -v expwin_i=$i '{print $1,$2-expwin_i,$3-expwin_i}' > 'HPC7_1.ER4_0.uniq.OnTad.joint.boundary.up'$i'.bed'
cat HPC7_1.ER4_0.uniq.OnTad.joint.boundary.bed | awk -F '\t' -v OFS='\t' -v expwin_i=$i '{print $1,$2+expwin_i,$3+expwin_i}' > 'HPC7_1.ER4_0.uniq.OnTad.joint.boundary.down'$i'.bed'
cat HPC7_1.ER4_1.uniq.OnTad.joint.boundary.bed | awk -F '\t' -v OFS='\t' -v expwin_i=$i '{print $1,$2-expwin_i,$3-expwin_i}' > 'HPC7_1.ER4_1.uniq.OnTad.joint.boundary.up'$i'.bed'
cat HPC7_1.ER4_1.uniq.OnTad.joint.boundary.bed | awk -F '\t' -v OFS='\t' -v expwin_i=$i '{print $1,$2+expwin_i,$3+expwin_i}' > 'HPC7_1.ER4_1.uniq.OnTad.joint.boundary.down'$i'.bed'
### up
bedtools intersect -a 'HPC7_1.ER4_0.uniq.OnTad.joint.boundary.up'$i'.bed' -b HPC7.13.ER4.13.mergedpk.bed -c > 'HPC7_1.ER4_0.HPC7_13.ER4_13.uniq.OnTad.joint.boundary.up'$i'.bed'
bedtools intersect -a 'HPC7_1.ER4_0.uniq.OnTad.joint.boundary.up'$i'.bed' -b HPC7.13.ER4.7.mergedpk.bed -c > 'HPC7_1.ER4_0.HPC7_13.ER4_7.uniq.OnTad.joint.boundary.up'$i'.bed'
bedtools intersect -a 'HPC7_1.ER4_1.uniq.OnTad.joint.boundary.up'$i'.bed' -b HPC7.13.ER4.13.mergedpk.bed -c > 'HPC7_1.ER4_1.HPC7_13.ER4_13.uniq.OnTad.joint.boundary.up'$i'.bed'
bedtools intersect -a 'HPC7_1.ER4_1.uniq.OnTad.joint.boundary.up'$i'.bed' -b HPC7.13.ER4.7.mergedpk.bed -c > 'HPC7_1.ER4_1.HPC7_13.ER4_7.uniq.OnTad.joint.boundary.up'$i'.bed'
### down
bedtools intersect -a 'HPC7_1.ER4_0.uniq.OnTad.joint.boundary.down'$i'.bed' -b HPC7.13.ER4.13.mergedpk.bed -c > 'HPC7_1.ER4_0.HPC7_13.ER4_13.uniq.OnTad.joint.boundary.down'$i'.bed'
bedtools intersect -a 'HPC7_1.ER4_0.uniq.OnTad.joint.boundary.down'$i'.bed' -b HPC7.13.ER4.7.mergedpk.bed -c > 'HPC7_1.ER4_0.HPC7_13.ER4_7.uniq.OnTad.joint.boundary.down'$i'.bed'
bedtools intersect -a 'HPC7_1.ER4_1.uniq.OnTad.joint.boundary.down'$i'.bed' -b HPC7.13.ER4.13.mergedpk.bed -c > 'HPC7_1.ER4_1.HPC7_13.ER4_13.uniq.OnTad.joint.boundary.down'$i'.bed'
bedtools intersect -a 'HPC7_1.ER4_1.uniq.OnTad.joint.boundary.down'$i'.bed' -b HPC7.13.ER4.7.mergedpk.bed -c > 'HPC7_1.ER4_1.HPC7_13.ER4_7.uniq.OnTad.joint.boundary.down'$i'.bed'
done



time Rscript plot_with_CI.R


HPC7_1.ER4_0.uniq.OnTad.joint.boundary.bed
14818
HPC7_1.ER4_1.uniq.OnTad.joint.boundary.bed
2607
HPC7.13.ER4.13.mergedpk.bed
4948
HPC7.13.ER4.7.mergedpk.bed
2905



