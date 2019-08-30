color_heatmap = function(color_matrix, outputname, format, border_color, w, h){
	format(outputname, width = w, height = h) ### output name
	par(mar=c(10,0.5,0.5,10)) ### set heatmap margins
	colbin_len = 10 ### column bin size
	rowbin_len = 10 ### row bin size
	### row reverse
	color_matrix = color_matrix[nrow(color_matrix):1,]
	### plot areas
	plot(c(0, dim(color_matrix)[2]*colbin_len), c(0, dim(color_matrix)[1]*rowbin_len), xaxt = "n", yaxt = "n", xaxs="i", yaxs="i", type = "n", xlab = "", ylab = "",main = "")
	### add color matrix colname as heatmap colname
	axis(1, c(1 : dim(color_matrix)[2])*colbin_len-0.5*colbin_len, colnames(color_matrix), las = 2, col.axis = "black", tick=FALSE)
	axis(4, c(1 : dim(color_matrix)[1])*colbin_len-0.5*colbin_len, rownames(color_matrix), las = 2, col.axis = "black", tick=FALSE)
	### use for loop to add rectangle with different color
	for (coln in c(1 : dim(color_matrix)[2])){ ### loop columns
		for (rown in c(1 : dim(color_matrix)[1])){ ### loop rows
			### add rectangle
			rect( (coln-1)*colbin_len, (rown-1)*rowbin_len, coln*colbin_len, rown*rowbin_len, col = color_matrix[rown, coln], border=border_color, lwd = 0 )
		}
	}
	dev.off()
}

mm10.10KB.chr19.CPneg.bed
m10.10KB.chr19.CPpos.bed
### read read counts
d_rc = read.table('/storage/home/g/gzx103/group/HiC/mouse/processed/10kb/G1E-ER4all3.merged.chr19.10kb.matrix', header=F)

color_mat = c()
#1:dim(d0)[1]
plot_region = c(4000:4500)
#plot_region = c(5000:6500)
d_rc_lim = d_rc
upperlim = 10
lowerlim = 0
d_rc_lim[d_rc_lim>upperlim] = upperlim
d_rc_lim[d_rc_lim<lowerlim] = lowerlim
lim_u = upperlim #max(d0_sub)
lim_l = lowerlim #min(d0_sub)
for (i in plot_region){
print(i)
color_row = c()
for (j in plot_region){
sig = d_rc_lim[i,j]
sig_p = (sig-lim_l)/(lim_u-lim_l)
r = ((0-255)*sig_p+255)/255
g= ((0-255)*sig_p+255)/255
b = ((0-255)*sig_p+255)/255
rgb_color = rgb(r,g,b)
color_row[j-min(plot_region)+1] = rgb_color
}
color_mat = rbind(color_mat, color_row)
}
colnames(color_mat) = rep('', length(plot_region))
rownames(color_mat) = rep('', length(plot_region))


### get rgb color table
color = read.table('colors.txt', header=F)

### get color matrix
color_mat_ideas = color_mat

type_cp = c('CPpos', 'CPneg')
color_cp = c('139,0,0', '16,78,139')

### get state
for ( si in 1:2){
print(si)
bed = read.table(paste('mm10.10KB.chr19.', type_cp[si], '.bed', sep=''), header=F)
print(dim(bed))
print(as.character(color_cp[si]))
state_rgb_code = as.numeric(unlist(strsplit(as.character(color_cp[si]), split=',')))
for (i in 1:dim(bed)[1]){
print(i)
#print(bed[i,])
for (j in 1:dim(bed)[1]){
if ((bed[i,4]-min(plot_region)+1>0) & (bed[j,4]-min(plot_region)+1>0) & (bed[i,4]-min(plot_region)+1<length(plot_region)) & (bed[j,4]-min(plot_region)+1<length(plot_region))){
sig = d_rc_lim[bed[i,4],bed[j,4]]
sig_p = (sig-lim_l)/(lim_u-lim_l)
r = ((state_rgb_code[1]-255)*sig_p+255)/255
g= ((state_rgb_code[2]-255)*sig_p+255)/255
b = ((state_rgb_code[3]-255)*sig_p+255)/255
rgb_color = rgb(r,g,b)
color_mat_ideas[bed[i,4]-min(plot_region)+1, bed[j,4]-min(plot_region)+1] = rgb_color
}
}
}
}


color_heatmap(color_mat_ideas, 'mm10.10KB.chr19.all.rc.CP.png', png, NA, 2000, 2000)





