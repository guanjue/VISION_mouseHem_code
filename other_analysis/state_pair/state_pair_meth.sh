#PBS -l nodes=1:ppn=8
#PBS -l walltime=20:00:00
#PBS -j oe
#PBS -A yzz2_e_g_sc_default

module load gcc/5.3.1
module load python/2.7.14-anaconda5.0.1
module load bedtools

cd /storage/home/gzx103/group/projects/vision/withMETH/check_order

cat /storage/home/g/gzx103/group/projects/vision/pknorm_2_16lim_ref1mo_sample1mo_0424/pknorm_2_16lim_ref1mo_0424_lesshet.state | awk -F ' ' -v OFS=' ' '{print $1,$2,$3,$4, $16,$9,$6,$12,$7,$18,$14, $25}' > noMeth.state


time Rscript get_state_pair.R
 




