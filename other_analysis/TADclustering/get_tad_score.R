chr_id_vec = c(1:19, 'X')
bed_tads_mat = c()

for (chr_id in chr_id_vec){
print(chr_id)
tadlevel2rgb = '../../TADlevel2rgb.txt'
input_tad_score_file = paste('OnTADraw_pen0.1_max200_meannorm_chr', chr_id, '_data1.tad', sep='')
input_tad_bed_file = paste('OnTADraw_pen0.1_max200_meannorm_chr', chr_id, '.bed', sep='')
ouput_file = paste('OnTADraw_pen0.1_max200_meannorm_chr', chr_id, '.withTADscore.bed', sep='')
### read tad score
d0 = read.table(input_tad_score_file, header=F)
d1 = d0[-1,]
d2 = d1[d1[,5]==0,6]
### read bed
bed = read.table(input_tad_bed_file, header=F)
### merge bed and tad score
bed_tads = cbind(bed, d2)
### get TAD level to rgb color
l2rgb = read.table(tadlevel2rgb, header=F)
### get TAD level from rgb
tadlevel = apply(bed_tads, 1, function(x) l2rgb[l2rgb[,2]==x[9],1])
### merge bed and tad score and tad level
bed_tads = cbind(bed, d2, tadlevel)
### write output
write.table(ouput_file, sep='\t', col.names=FALSE, row.names=FALSE, quote=FALSE)
###
bed_tads_mat = rbind(bed_tads_mat, bed_tads)
}

write.table(bed_tads_mat, 'bed_tads_bed_withscore.bed', sep='\t', col.names=FALSE, row.names=FALSE, quote=FALSE)
