cat cmp9_eryad12.bed | awk -F '\t' -v OFS='\t' '{if (int(($2+$3)/2-5000)>0) print $1, int(($2+$3)/2-5000), int(($2+$3)/2+5000); else print $1, 0, int(($2+$3)/2+5000)}' > cmp9_eryad12.5kbupdown.bed
cat cmp9_eryad12.bed | awk -F '\t' -v OFS='\t' '{if (int(($2+$3)/2-5000)>0) print $1, int(($2+$3)/2-5000), int(($2+$3)/2+5000)}' > cmp9_eryad12.5kbupdown.bed
bedtools getfasta -fi ~/group/genome/mm10/mm10_no_alt_analysis_set_ENCODE.fasta -fo cmp9_eryad12.5kbupdown.bed.fa -bed cmp9_eryad12.5kbupdown.bed

cat cmp9_eryad3.bed | awk -F '\t' -v OFS='\t' '{if (int(($2+$3)/2-5000)>0) print $1, int(($2+$3)/2-5000), int(($2+$3)/2+5000); else print $1, 0, int(($2+$3)/2+5000)}' > cmp9_eryad3.5kbupdown.bed
cat cmp9_eryad3.bed | awk -F '\t' -v OFS='\t' '{if (int(($2+$3)/2-5000)>0) print $1, int(($2+$3)/2-5000), int(($2+$3)/2+5000)}' > cmp9_eryad3.5kbupdown.bed
bedtools getfasta -fi ~/group/genome/mm10/mm10_no_alt_analysis_set_ENCODE.fasta -fo cmp9_eryad3.5kbupdown.bed.fa -bed cmp9_eryad3.5kbupdown.bed


time ~/group/software/meme/bin/fimo --text --thresh 1 cmp9_ery12.sequnwinder.motif1.meme cmp9_eryad12.5kbupdown.bed.fa > cmp9_eryad12.5kbupdown.gata.fimo.txt
time ~/group/software/meme/bin/fimo --text --thresh 1 cmp9_ery12.sequnwinder.motif1.meme cmp9_eryad3.5kbupdown.bed.fa > cmp9_eryad3.5kbupdown.gata.fimo.txt
time ~/group/software/meme/bin/fimo --text --thresh 1 cmp9_ery3.sequnwinder.motif1.meme cmp9_eryad12.5kbupdown.bed.fa > cmp9_eryad12.5kbupdown.etf.fimo.txt
time ~/group/software/meme/bin/fimo --text --thresh 1 cmp9_ery3.sequnwinder.motif1.meme cmp9_eryad3.5kbupdown.bed.fa > cmp9_eryad3.5kbupdown.etf.fimo.txt


time python find_best_match_fimo.py -i cmp9_eryad12.5kbupdown.gata.fimo.txt -o cmp9_eryad12.5kbupdown.gata.fimo.maxpos.txt
time python find_best_match_fimo.py -i cmp9_eryad3.5kbupdown.gata.fimo.txt -o cmp9_eryad3.5kbupdown.gata.fimo.maxpos.txt
time python find_best_match_fimo.py -i cmp9_eryad12.5kbupdown.etf.fimo.txt -o cmp9_eryad12.5kbupdown.etf.fimo.maxpos.txt
time python find_best_match_fimo.py -i cmp9_eryad3.5kbupdown.etf.fimo.txt -o cmp9_eryad3.5kbupdown.etf.fimo.maxpos.txt

time Rscript get_motif_max_pos_heatmap.R cmp9_eryad12.5kbupdown.gata.fimo.maxpos.txt cmp9_eryad12.5kbupdown.gata.fimo.maxpos.png 10000
time Rscript get_motif_max_pos_heatmap.R cmp9_eryad3.5kbupdown.gata.fimo.maxpos.txt cmp9_eryad3.5kbupdown.gata.fimo.maxpos.png 10000
time Rscript get_motif_max_pos_heatmap.R cmp9_eryad12.5kbupdown.etf.fimo.maxpos.txt cmp9_eryad12.5kbupdown.etf.fimo.maxpos.png 10000
time Rscript get_motif_max_pos_heatmap.R cmp9_eryad3.5kbupdown.etf.fimo.maxpos.txt cmp9_eryad3.5kbupdown.etf.fimo.maxpos.png 10000






.libPaths("/storage/home/gzx103/R/x86_64-redhat-linux-gnu-library")

> .libPaths()
[1] "/opt/aci/sw/r/3.4/lib64/library_gcc-5.3.1"
[2] "/opt/aci/sw/r/3.4/share/library"          
[3] "/usr/lib64/R/library"                     
[4] "/usr/share/R/library"

