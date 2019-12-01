### get lm coefficients
folder_of_Rdat='/gpfs/group/yzz2/default/legacy/group/projects/vision/rna/'
time Rscript regress_coefficient_2kb.with0.R $folder_of_Rdat

### get cCRE state proportion
cd /storage/home/gzx103/group/projects/vision/pknorm_2_16lim_ref1mo_sample1mo_0424/get_eRP_tracks
### get state counts in cCRE in each ct
qsub PBS.splitbed.sh

### cbind state count
time bash split.sh

### get 8 track (4 gene groups, proximal & distal) for each cCRE 
qsub PBS.getbw.sh
