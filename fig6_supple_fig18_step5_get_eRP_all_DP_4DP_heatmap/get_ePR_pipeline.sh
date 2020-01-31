cd ~/group/projects/vision/rna

### get eRPs raw coefficient
time Rscript regress_eRP_2kb.R

### get eRPs raw coefficient -> split D & P
time Rscript regress_eRP_2kb_DP.R

### get eRPs split D & P -> split 4 groups
time Rscript regress_eRP_2kb_split4group.R

### plot eRPs
time Rscript plot_eRPs.R

