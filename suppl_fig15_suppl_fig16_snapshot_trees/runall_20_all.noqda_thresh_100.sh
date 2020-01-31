#PBS -l nodes=1:ppn=8
#PBS -l walltime=10:00:00
#PBS -j oe
#PBS -A yzz2_e_g_sc_default
#PBS -l pmem=16gb

module load gcc/5.3.1 
module load python/2.7.14-anaconda5.0.1
module load bedtools

cd /storage/home/gzx103/scratch/vision/5end/index_set_20cell_pknorm_qda


##################################
script_folder='/storage/home/gzx103/group/software/snapshot/bin/'
input_folder='/storage/home/gzx103/scratch/vision/5end/index_set_20cell_pknorm_qda/'
output_folder='/storage/home/gzx103/scratch/vision/5end/index_set_20cell_pknorm_qda/snapshot18_reproduce_0_16lim_all_0/'
qda_num=0
#mkdir atac_sig
#cp /storage/home/gzx103/group/projects/vision/pknorm_16lim_ref1mo_sample1mo_correct/*.atac.pkn16.txt atac_sig/
#cp /storage/home/gzx103/scratch/vision/5end/pknorm_16lim/pcor_100lim_mean/pknorm_16lim_ref1mo_sample1mo_correct/*.atac.pkn16.txt atac_sig/
#ls atac_sig/*.atac.pkn16.txt > atac_sig_list.txt
#for file in $(cat atac_sig_list.txt)
#do
#        paste 200bp.win.bed $file > $file'.bed'
#done

#mkdir ideas_label
#cp /storage/home/gzx103/group/projects/vision/pknorm_2_16lim_ref1mo_sample1mo_0424/pknorm_2_16lim_ref1mo_0424_lesshet.state ideas_label/used.state
#head -1 ideas_label/used.state > ideas_label/used.modified.state
#tail -n+2 ideas_label/used.state >> ideas_label/used.modified.state

#cd ideas_label/
#time python /storage/home/gzx103/group/software/CD_viewer/bin/ideas_matrix2ideas_bed.py -i used.state -a 2 -b 5 -c 5 -d 24 -o 20ct
#cd ..
#ls ideas_label/20ct.*.bed > ideas_list.txt
#for file in $(cat atac_sig_list.txt)
#do
#        cat $file | awk -F '\t' -v OFS='\t' '{print $1,$2,$3,".",$4}' > $file'.tmp.txt'
#        mv $file'.tmp.txt' $file
#done

### run snapshot (CORE!!!)
echo 'run snapshot :o'
time python $script_folder'snapshot_v.0.3.py' -p peak_list.txt -n atac_20cell -t 100 -s signal_list.txt -l F -z F -x 0.01 -f ideas_list.txt -m mostfreq -c ideas_range_color.txt -e cd_tree.txt -i $input_folder -o $output_folder -b $script_folder -q $qda_num
#time python 'snapshot_v.0.1.py' -p peak_list.txt -n atac_20cell -t 100 -s signal_list.txt -l F -z F -x 0.01 -f ideas_list.txt -m mostfreq -c ideas_range_color.txt -e cd_tree.txt -i $input_folder -o $output_folder -b $script_folder
echo 'complete :)'

cd snapshot20_reproduce
#paste atac_20cell.index.matrix.txt atac_20cell.function.matrix.txt | awk -F '\t' -v OFS='\t' '{print $1,$2,$3, $5"_"$6"_"$7"_"$8"_"$9"_"$10"_"$11"_"$12"_"$13"_"$14"_"$15"_"$16"_"$17"_"$18"_"$19"_"$20"_"$21"_"$22, $27,$28,$29,$30,$31,$32,$33,$34,$35,$36,$37,$38,$39,$40,$41,$42,$43,$44,$45,$46}' > atac_20cell.fun.txt
#paste atac_20cell.index.matrix.no0.txt atac_20cell.function.matrix.no0.txt | awk -F '\t' -v OFS='\t' '{print $1,$2,$3, $5"_"$6"_"$7"_"$8"_"$9"_"$10"_"$11"_"$12"_"$13"_"$14"_"$15"_"$16"_"$17"_"$18"_"$19"_"$20"_"$21"_"$22, $27,$28,$29,$30,$31,$32,$33,$34,$35,$36,$37,$38,$39,$40,$41,$42,$43,$44,$45,$46}' > atac_20cell.fun.no0.txt
cd ..

