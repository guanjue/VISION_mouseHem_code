### get lm coefficients
folder_of_Rdat='/gpfs/group/yzz2/default/legacy/group/projects/vision/rna/'
time Rscript regress_coefficient_2kb.with0.R $folder_of_Rdat

### get state counts in cCRE in each ct
qsub PBS.splitbed.sh

### cbind state count
time bash split.sh

### get 8 track (4 gene groups, proximal & distal) for each cCRE 
qsub PBS.getbw.sh
