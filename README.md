# VISION mouse pipeline: from  Reads Count (RC) to PKnorm normalized -log10(p-value)

#### Advanced sequencing technologies have generated a plethora of data for many chromatin marks in multiple tissues and cell types, yet there is lack of a generalized tool for optimal utility of those data. A major challenge is to quantitatively model the epigenetic dynamics across both the genome and many cell types for understanding their impacts on differential gene regulation and disease. We introduce IDEAS, an integrative and discriminative epigenome annotation system, for jointly characterizing epigenetic landscapes in many cell types and detecting differential regulatory regions. A key distinction between our method and existing state-of-the-art algorithms is that IDEAS integrates epigenomes of many cell types simultaneously in a way that preserves the position-dependent and cell type-specific information at fine scales, thereby greatly improving segmentation accuracy and producing comparable annotations across cell types. 
#### (Zhang, Yu, Lin An, Feng Yue, and Ross C. Hardison. "Jointly characterizing epigenetic dynamics across multiple human cell types." Nucleic acids research 44, no. 14 (2016): 6721-6731.)


<img src="https://github.com/guanjue/IDEAS_2018/blob/master/example_figures/" width="800"/>

##### Figure 1. Illustration of IDEAS model. The IDEAS method borrows locus specific information across cell types to improve accuracy, and simultaneously accounts for local cell type relationships for inferring cell type-specific activities. In particular, to infer epigenetic state at a given locus in a target cell type, IDEAS uses the currently inferred states in other cell types at the same locus, but only those cell types showing similar local epigenetic landscapes with the target cell type, as priors to improve inference. The local window (dashed box) is dynamically determined by Markov chains, and cell types are clustered within the local window for their relationships with the target cell. The entire process is iterative with all the unknowns (epigenetic states and local cell type clustering) updated until convergence. The final segmentation is then colored using an automatic coloring script for visualization in browser.




## Install IDEAS
#### Clone the github repository 
```
git clone https://github.com/guanjue/vision_mouse.git
```


## Input data
##### The parameter file for IDEAS. 
```
run_IDEAS.parafile
```
##### Usually, user only needs to change the following parameters in the parameter file:
```
script_dir
working_dir
input_file_list
input_dir
overall_upper
overall_lower
```


##### The bin file for IDEAS: each column is separated by whitespace
###### 1st column: target dataset; 
###### 2nd column: the no antibody control file for the target dataset; 
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
## Run IDEAS
##### (1) copy the 'run_pipeline.sh' into the working directory
```
cp ~/group/software/IDEAS/IDEAS_2018/run_pipeline.sh working_dir/

```
##### (2) change the following parameters in the 'run_IDEAS.sh' file:
###### script_dir='absolute path to the IDEAS_2018 dir'
###### working_dir='absolute path to the working directory'
###### input_file_list='the input and the correponding no antibody control file list'
###### input_dir='the input file's folder'
###### overall_upper='the upper limit of the output file'
###### overall_lower='the lower limit of the output file'

```
>>> head -100 run_IDEAS.sh 
###### set parameters
script_dir=/storage/home/gzx103/group/software/vision_mouse/src/
working_dir=/storage/home/gzx103/scratch/vision/test_pipeline/
input_dir=/storage/home/gzx103/scratch/vision/test_pipeline/input_5end_rc/
input_file_list=info_table_all.rc2nbp.txt
overall_upper=16
overall_lower=2


time bash $script_dir'overall_pipeline.sh' $script_dir $working_dir $input_dir $input_file_list $overall_upper $overall_lower
```

##### (4) use 'run_IDEAS.sh' script to run IDEAS
```
time bash run_pipeline.sh
```



## Output results for test data
### All output files will be saved to the following directories inside the working directory:
```
nbp/
pknorm_info/
pknorm_sig/
ref_info/
pknorm_ref_sig/
fisherp/
list_files/
pknorm_2_16_sig/
```

## The heatmap for IDEAS epigenetic state
<img src="https://github.com/guanjue/IDEAS_2018/blob/master/example_figures/f3_vision_result.png" width="800"/>

##### Figure 3. The epigenetic state inferred by IDEAS. Each row represents one epigenetic state. Each column represents one epigenetic mark. The color density represent the average signal of all genome loci with the corresponding epigenetic state. The dark blue represents high average signal. The white represent the low average signal. 

## The genome browser track (bigbed format) for IDEAS epigenetic state will be saved in the subdirectory (named as Tracks/) in the output directory: 
```
track_dir=/storage/home/gzx103/group/software/IDEAS/IDEAS_2018/test_data/run_IDEAS_result/Tracks/
```

## References

##### Zhang, Yu, Lin An, Feng Yue, and Ross C. Hardison. "Jointly characterizing epigenetic dynamics across multiple human cell types." Nucleic acids research 44, no. 14 (2016): 6721-6731.
##### Zhang, Yu, and Ross C. Hardison. "Accurate and reproducible functional maps in 127 human cell types via 2D genome segmentation." Nucleic acids research 45, no. 17 (2017): 9823-9836.


