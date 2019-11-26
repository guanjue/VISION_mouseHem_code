#PBS -l nodes=1:ppn=8
#PBS -l walltime=10:00:00
#PBS -j oe
#PBS -A yzz2_e_g_sc_default
#PBS -l pmem=16gb

module load gcc
module load python/2.7.14-anaconda5.0.1


cd /storage/home/g/gzx103/scratch/vision/5end/pknorm_16lim/pcor_100lim_mean/pknorm_100lim_ref1mo_sample1mo

time Rscript plot_correlation_heatmap_withcolorbar.R pknorm_100lim_ref1mo_sample1mo_list.txt red white black average pearson pr_cor_heat_a.pdf

time Rscript plot_correlation_heatmap_withcolorbar.R pknorm_100lim_ref1mo_sample1mo_list.txt red white black average spearman sp_cor_heat_a.pdf

