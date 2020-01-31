time Rscript plot_all_fun_bar.R atac_pk_no0.atac.ideas.txt atac_pk_no0.atac.index.txt signal_list_atac.txt ideas_range_color.txt atac_18cell.allfunbar

time Rscript plot_all_fun_bar_noneedpk.R atac_pk_no0.atac.ideas.txt atac_pk_no0.atac.index.txt signal_list_atac.txt ideas_range_color.txt atac_18cell.allfunbar

time Rscript plot_all_fun_bar_s3normpk.R atac_pk_no0.atac.ideas.txt atac_pk_no0.atac.signal.txt signal_list_atac.txt ideas_range_color.txt atac_18cell.allfunbar 3

### get pk num & ideas state num
for i in {4..7}
do
echo $i
time Rscript plot_all_fun_bar_s3normpk.R atac_pk_no0.atac.ideas.txt atac_pk_no0.atac.signal.txt signal_list_atac.txt ideas_range_color.txt atac_18cell.allfunbar $i
done

### get hist peak enrichment
time Rscript time Rscript atac_vs_histone.R


### 0804 fun bar only 9 & 13
time Rscript plot_all_fun_bar.only9_13.R atac_pk_no0.atac.ideas.txt atac_pk_no0.atac.index.txt signal_list_atac.txt ideas_range_color.txt atac_18cell.allfunbar
### 0804 fun bar not 0 & 9 & 13
time Rscript plot_all_fun_bar.not0_9_13.R atac_pk_no0.atac.ideas.txt atac_pk_no0.atac.index.txt signal_list_atac.txt ideas_range_color.txt atac_18cell.allfunbar
