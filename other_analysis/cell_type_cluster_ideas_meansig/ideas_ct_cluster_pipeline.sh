script_dir='/storage/home/gzx103/group/software/vision_mouse/other_analysis/cell_type_cluster_ideas_meansig/'

######
### get ideas mean signal heatmap & mean signal matrix
time Rscript $script_dir'plot_ideas_heatmap.R'


######
### get mv distance of each bin based on IDEAS state mean signal
time python $script_dir'ct_dist_ideas_state_mean_signal.py' -i pknorm_2_16lim_ref1mo_0424_lesshet.state -s ideas_state_mean_signal.txt -o IDEAS_state_meansignal_Euclidean_dist_matrix.txt


######
### plot mv Euclidean distance heatmap 
time Rscript $script_dir'plot_ct_cluster_heatmap.R'
