# VISION mouse pipeline: from  Reads Count (RC) to S3norm normalized -log10(p-value)

#### Here, we developed a data preprocessing pipeline for normalizing the sequencing depth and the signal-to-noise ratio across multiple cell types in the ValIdated Systematic IntegratiON of hematopoietic epigenomes (VISION) data consortium project. The advantage of this pipeline is that it can normalize the signal of the peak regions without expand the signal of the background regions across multiple dataset that were generated by different labs.

#### VISION project dataset portal: usevision.org


<img src="https://github.com/guanjue/vision_mouse/blob/master/example_figures/rc2pknorm.png" width="800"/>

##### Figure 1. Illustration of the signal track before and after data preprocessing pipeline. In this pipeline, we started with the reads count of 5'-end in each 200-bp bin and end with S3norm normalized -log10(p-value) based on a dynamic negative binomial (NB) model. In the pipeline, we first counted the number of reads in each of the 200-bp bin. Then, we used a two rounds negative binomial model to convert the reads count to the -log10(p-value). For replicates p-value of the same data, we used the Fisher's method to merge them. In the third step, we used the S3norm normalization to normalize both the sequencing depth and the the signal-to-noise ratio of these fisher-merged-p-values. In the third step, we first selected the dataset with highest Fraction of signal in peaks (FRiP score) as the reference dataset for each mark. These reference dataset were normalized by a modified S3norm method (matching the mean signal of all peaks in each two datasets, instead of only matching mean signal of the common peak regions). Then, we used the cross mark normalized reference dataset as the reference dataset to normalize each dataset with the same mark. To avoid outlier effects, we limit our data by set a upper limit equals to 16 (Left) Before the data preprocessing, the signal of the peak regions are not comparable across all dataset. (Middle) We first convert the reads count to -log10(p-value) based on a dynamic negative binomial background model. (Right) After the S3norm the signal of the peak regions and the background regions become more comparable across all datasets.




## Install VISION mouse pipeline
#### Clone the github repository 
```
git clone https://github.com/guanjue/vision_mouse.git
```
#### Install dependency: go to vision_mouse folder and run the following command
```
time bash INSTALL.sh
```



##### The input file list for VISION mouse pipeline: 
###### 1st column: the signal of each bin;
###### all of the input file should be saved in the input folder (input_dir)
```
info_table_all.rc2nbp.txt
>>> head -1000000 B_SPL.h3k27acrep.100035.bamtobed5endintersect.signal | tail -20
0
0
0
1
1
0
0
0
2
0
0
0
0
4
```

##### The input filename list for VISION mouse pipeline: each column is separated by tab
###### 1st column: target dataset; 
###### 2nd column: the no antibody control file for the target dataset; Here, we used the same merge input file base on the 21 input files from 11 cell types. For each input file, it is first normalized to ER4.137.input.signal input dataset based the ratio of total reads count. Then, the same merge input file is the mean signal of the 21 normalized the input files. 
###### For the atac-seq without input signal. We used a input file with all bins equal to one as the input signal.
```
info_table_all.rc2nbp.txt
>>> head info_table_all.rc2nbp.txt
B_SPL.h3k27acrep.100035.bamtobed5endintersect.signal	merged_normed_input.rounding.txt
CFU_E_ad.h3k27acrep.100030.bamtobed5endintersect.signal	merged_normed_input.rounding.txt
CLP.h3k27acrep.100039.bamtobed5endintersect.signal	merged_normed_input.rounding.txt
CMP.h3k27acrep.100027.bamtobed5endintersect.signal	merged_normed_input.rounding.txt
ER4.h3k27acrep.538.bamtobed5endintersect.signal	merged_normed_input.rounding.txt
ER4.h3k27acrep.539.bamtobed5endintersect.signal	merged_normed_input.rounding.txt
......
```


## Run VISION mouse pipeline
##### (1) copy the 'run_pipeline.sh' in the VISION mouse pipeline into the working directory
```
cp ~/group/software/vision_mouse/run_pipeline.sh working_dir/
```
##### (2) change the following parameters in the 'run_IDEAS.sh' file:
###### script_dir='absolute path to the IDEAS_2018 dir'
###### working_dir='absolute path to the working directory'
###### input_file_list='the input and the correponding no antibody control file list'
###### input_dir='the input file's folder'
###### overall_upper='the upper limit of the output file'
###### overall_lower='the lower limit of the output file'
###### select_method='method used to select reference dataset (frip/snr)'
###### user_given_global_ref='user given global reference dataset (if empty, pipeline will user the dataset with the highest frip/snr score dataset)'

```
>>> head -100 run_IDEAS.sh 
###### set parameters
script_dir=/storage/home/gzx103/group/software/vision_mouse/src/
working_dir=/storage/home/gzx103/scratch/vision/test_pipeline/
input_dir=/storage/home/gzx103/scratch/vision/test_pipeline/input_5end_rc/
input_file_list=info_table_all.rc2nbp.txt
overall_upper=16
overall_lower=2
select_method=frip
user_given_global_ref=ERY_fl.h3k4me3rep.fisher_p.txt

time bash $script_dir'overall_pipeline.sh' $script_dir $working_dir $input_dir $input_file_list $overall_upper $overall_lower $select_method $user_given_global_ref
```

##### (3) use 'run_IDEAS.sh' script to run VISION mouse pipeline
```
time bash run_pipeline.sh
```



## Output results for test data
### All output files will be saved to the following directories inside the working directory:
```
nbp/: The p-value of each sample based on the NB background model.
pknorm_info/: The scatterplot and the parameters used in S3norm normalization for each dataset
pknorm_sig/: The signal files of the S3norm normalized data
ref_info/: The scatterplot and the parameters used in S3norm normalization for reference datasets across all marks
pknorm_ref_sig/: The signal files of the S3norm normalized reference data
fisherp/: The merged p-value from each sample's NB p-value by the Fisher's method
list_files/: All of the list files used in the pipeline
pknorm_2_16_sig/: The signal files of the S3norm normalized data with the upper & lower bound limitation
```


## References

##### PKnorm & vision paper
Xiang, Guanjue, et al. "S3norm: simultaneous normalization of sequencing depth and signal-to-noise ratio in epigenomic data." bioRxiv (2019): 506634.

Xiang, Guanjue, et al. "An integrative view of the regulatory and transcriptional landscapes in mouse hematopoiesis." bioRxiv (2020): 731729.


# vision_mouse_mm10
