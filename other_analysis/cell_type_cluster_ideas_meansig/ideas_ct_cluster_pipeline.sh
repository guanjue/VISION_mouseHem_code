script_dir=

######
### get ideas mean signal heatmap & mean signal matrix
time Rscript plot_ideas_heatmap.R


######
### get mv distance of each bin based on IDEAS state mean signal
time python ct_dist_ideas_state_mean_signal.py -i pknorm_2_16lim_ref1mo_0424_lesshet.state -s ideas_state_mean_signal.txt -o IDEAS_state_meansignal_Euclidean_dist_matrix.txt

