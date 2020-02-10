### get nb-p-value
cd /storage/home/gzx103/scratch/vision/5end/
cp h3k27ac/h3k27ac_nbp_bgadj.sh /storage/home/gzx103/scratch/vision/used_scripts/nbp_bgadj/
cp h3k27ac/info_table_histone_input.txt /storage/home/gzx103/scratch/vision/used_scripts/nbp_bgadj/info_table_histone_input.h3k27ac.txt
cp h3k27me3/h3k27me3_nbp_bgadj.sh /storage/home/gzx103/scratch/vision/used_scripts/nbp_bgadj/
cp h3k27me3/info_table_histone_input.txt /storage/home/gzx103/scratch/vision/used_scripts/nbp_bgadj/info_table_histone_input.h3k27me3.txt
cp h3k36me3/h3k36me3_nbp_bgadj.sh /storage/home/gzx103/scratch/vision/used_scripts/nbp_bgadj/
cp h3k36me3/info_table_histone_input.txt /storage/home/gzx103/scratch/vision/used_scripts/nbp_bgadj/info_table_histone_input.h3k36me3.txt
cp h3k9me3/h3k9me3_nbp_bgadj.sh /storage/home/gzx103/scratch/vision/used_scripts/nbp_bgadj/
cp h3k9me3/info_table_histone_input.txt /storage/home/gzx103/scratch/vision/used_scripts/nbp_bgadj/info_table_histone_input.h3k9me3.txt
cp h3k4me3/h3k4me3_nbp_bgadj.sh /storage/home/gzx103/scratch/vision/used_scripts/nbp_bgadj/
cp h3k4me3/info_table_histone_input.txt /storage/home/gzx103/scratch/vision/used_scripts/nbp_bgadj/info_table_histone_input.h3k4me3.txt
cp h3k4me1/h3k4me1_nbp_bgadj.sh /storage/home/gzx103/scratch/vision/used_scripts/nbp_bgadj/
cp h3k4me1/info_table_histone_input.txt /storage/home/gzx103/scratch/vision/used_scripts/nbp_bgadj/info_table_histone_input.h3k4me1.txt
cp ctcf/ctcf_nbp_bgadj.sh /storage/home/gzx103/scratch/vision/used_scripts/nbp_bgadj/
cp ctcf/info_table_tf_input.txt /storage/home/gzx103/scratch/vision/used_scripts/nbp_bgadj/info_table_tf_input.txt
cp atac/atac_nbp_bgadj.sh /storage/home/gzx103/scratch/vision/used_scripts/nbp_bgadj/
cp atac/atac_list.txt /storage/home/gzx103/scratch/vision/used_scripts/nbp_bgadj/info_table_atact_input.txt
### input data
cp h3k27ac/merged_normed_input.rounding.txt /storage/home/gzx103/scratch/vision/used_scripts/nbp_bgadj/
cp ~/group/projects/vision/merged_input/ones_input.ataconly.rounding.txt /storage/home/gzx103/scratch/vision/used_scripts/

### get fisher-p
cd /storage/home/gzx103/scratch/vision/5end/fisher_p_100lim/
cp get_bw.sh /storage/home/gzx103/scratch/vision/used_scripts/fisherp/

### get pknorm
# ref
cd /storage/home/gzx103/scratch/vision/5end/pknorm_16lim/ref_r_log_m
cp pknorm_ref_r_refmean.sh /storage/home/gzx103/scratch/vision/used_scripts/pknorm/refnorm/
cp raw_sig_list.txt /storage/home/gzx103/scratch/vision/used_scripts/pknorm/refnorm/
# ct
cd /storage/home/gzx103/scratch/vision/5end/pknorm_16lim/
cp *_r_log_mm/pknorm_*_r.sh /storage/home/gzx103/scratch/vision/used_scripts/pknorm/

### get input_norm
cd /storage/home/gzx103/group/projects/vision/
cp merged_input/total_mean_norm_input.sh /storage/home/gzx103/scratch/vision/used_scripts/input_norm/

