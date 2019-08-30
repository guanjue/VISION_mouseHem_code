cd /storage/home/gzx103/group/projects/vision/raw_5end

qsub get_bedgraph.sh

paste ../bin.tab.mm10.200bp.noblacklist.bed ../merged_input/merged_normed_input_mean.rounding.txt | awk -F '\t' -v OFS='\t' '{print $1,$2,$3,$5}' | sort -k1,1 -k2,2n > merged_normed_input_mean.rouding.bedgraph



