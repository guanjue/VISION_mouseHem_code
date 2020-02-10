### step 1
### select gene-cCRE pair
cd vision_rna_tss2k_ccreunit_folder_used0509

### step 2
### extract gene cCRE from Rdat
cd /storage/home/gzx103/group/projects/vision/rna/
bash extract.gene_ccRE.sh
head vision_rna.gene_ccRE.selected.all.reprod_count.WithName1209.txt 

### step 3
### get matrix with chr tss gene_name ccRE_start ccRE_end correlation selection_num gene_group cCRE_states
cd /storage/home/gzx103/group/projects/vision/rna/gene_ccRE_200bp_to_ccRE/
bash gene_ccRE_200bp_to_ccRE.sh
head all_converted.sort.g.cCRE_state.txt

### step 4
### get 8 eRP scores track for each cCRE in each cell type
cd /storage/home/gzx103/group/projects/vision/pknorm_2_16lim_ref1mo_sample1mo_0424/get_eRP_tracks
bash PBS.splitbed.sh
bash split.sh
bash PBS.getbw.sh

### step 7
### add 8 eRP scores track for each cCRE in each cell type
cd ~/group/projects/vision/rna/ccRE_eRPs_Rdat_version/
bash get_ccRE_eRP_for_megatable.sh
head mouse_allChr_wTAD_022519.eRPs.txt



