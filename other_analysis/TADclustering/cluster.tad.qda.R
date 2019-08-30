library(pheatmap)
library(mclust)
library(MASS)
library(mixtools)

stateColor<-function(statemean, markcolor=NULL)
{	
	if(length(markcolor)==0)
	{	markcolor=rep("",dim(statemean)[2]);
		markcolor[order(apply(statemean,2,sd),decreasing=T)]=hsv((1:dim(statemean)[2]-1)/dim(statemean)[2],1,1)
		markcolor=t(col2rgb(markcolor));
	}
	rg=apply(statemean,1,range);
	mm=NULL;
	for(i in 1:dim(statemean)[1])
	{	mm=rbind(mm,(statemean[i,]-rg[1,i]+1e-10)/(rg[2,i]-rg[1,i]+1e-10));
	}
	mm = mm^5; 
	if(dim(mm)[2]>1) mm = mm / (apply(mm, 1, sum)+1e-10);
	#print(mm)
	mycol=mm%*%markcolor;
	s=apply(statemean,1,max);
	s=(s-min(s))/(max(s)-min(s)+1e-10);
#s=s^0.5;
mycol=round(255-(255-mycol)*s/0.8);
mycol[mycol<0]=0;
rt=paste(mycol[,1],mycol[,2],mycol[,3],sep=",");
h=t(apply(mycol,1,function(x){rgb2hsv(x[1],x[2],x[3])}));
h=apply(h,1,function(x){hsv(x[1],x[2],x[3])});
rt=cbind(rt,h);
return(rt);

	h=t(apply(mycol,1,function(x){rgb2hsv(x[1],x[2],x[3])}));
	h[,2]=h[,2]*s;
	#h[,3]=1-(1-h[,3])*s^0.5;
	h=apply(h,1,function(x){hsv(x[1],x[2],x[3])});
	rt=cbind(apply(t(col2rgb(h)),1,function(x){paste(x,collapse=",")}),h);
	
	return(rt);
}

color_heatmap = function(color_matrix, outputname, format, border_color, w, h){
	format(outputname, width = w, height = h) ### output name
	par(mar=c(12,0.5,0.5,12)) ### set heatmap margins
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
			rect( (coln-1)*colbin_len, (rown-1)*rowbin_len, coln*colbin_len, rown*rowbin_len, col = color_matrix[rown, coln], border=border_color, lwd = 0.1 )
		}
	}
	dev.off()
}

set.seed(2018)
### read signal matrix
d1_raw0 = c()
for (i in 2:13){
#for (i in 1:27){
	#filename = paste('IDEAS_input_merged_tad_and_bd/IDEAS_state_p.', i, '.txt', sep='')
	filename = paste('IDEAS_input_merged_tadbdonly/IDEAS_state_p.', i, '.txt', sep='')
	sig = scan(filename)
	d1_raw0 = cbind(d1_raw0, sig)
}

### get index
d1_raw = (d1_raw0)
#d1_raw_binary = d1_raw
#d1_raw_binary[d1_raw_binary<=0] = 0
#d1_raw_binary[d1_raw_binary>0] = 2
#d1_raw_binary[d1_raw_binary>=1 & d1_raw_binary<2] = 1

lim = c()
pdf('hist_sig.pdf', width=21, height=21)
par(mfrow=c(4,3))
for (i in 1:dim(d1_raw)[2]){
hist(log2(d1_raw[,i]), breaks=50)
}
dev.off()

#lim = c()
#pdf('hist_sig.gmm.pdf', width=21, height=14)
#par(mfrow=c(2,3))
#for (i in 1){
#mixmdl = normalmixEM((d1_raw[,i]), k=3)
#plot(mixmdl,which=2)
#lines(density((d1_raw[,i])), lty=2, lwd=2)
#lim = rbind(lim, c(max(d1_raw[mixmdl$posterior[,1]>0.5,i]), max(d1_raw[mixmdl$posterior[,1]<=0.5,i])))
#}
#dev.off()


#pdf('hist_sig.all.pdf', width=21, height=14)
#par(mfrow=c(1,1))
#mixmdl = normalmixEM(log2(d1_raw), k=3)
#plot(mixmdl,which=2)
#lines(density(log2(d1_raw[,i])), lty=2, lwd=2)
#dev.off()

#lim_1 = apply(lim, 1, min)

d1_raw_binary = d1_raw
for ( i in 1:dim(d1_raw)[2]){
	d1_raw_z_i = (d1_raw[,i]-mean(d1_raw[,i])) / sd(d1_raw[,i])
	d1_raw_zp_i = 1 - pnorm(d1_raw_z_i) 
	d1_raw_fdrzp_i = p.adjust(d1_raw_zp_i)
	d1_raw_binary[d1_raw_fdrzp_i>1e-3,i] = 0
	d1_raw_binary[d1_raw_fdrzp_i<=1e-3,i] = 1
}

### get ISs
d1_raw_binary_index = apply(d1_raw_binary, 1, function(x) paste(x, collapse = '_'))
count_threshold = 10
### get IS threshold
index_count = table(d1_raw_binary_index)
if (length(index_count)>100){
top = 0.95
mean_95 = mean(index_count[index_count<=quantile(index_count, top)])
var_95 = var(index_count[index_count<=quantile(index_count, top)])
sd_95 = sd(index_count[index_count<=quantile(index_count, top)])
size = (mean_95^2)/(var_95 - mean_95)
prob=mean_95/var_95
if (prob < 0.001){
	prob = 0.001
} else if (prob > 0.999){
	prob = 0.999
}
size = prob*mean_95 / (1-prob)
### get NB p-val
pvec = pnbinom(index_count, size=size, prob=prob, lower.tail = FALSE)
print(summary(as.matrix(pvec)))
### FDR pval
padjvec = p.adjust(pvec, method='bonferroni')
print(summary(as.matrix(padjvec)))
### get NB count thresh
counts_pfdr = cbind(index_count, padjvec)
if (sum(padjvec<0.01)>0){
	print('sum(padjvec<0.01)>0')
	NB_count_thresh = min(counts_pfdr[padjvec<0.01,1])
} else {
	print('user provide')
	NB_count_thresh = count_threshold
	print('NB model fail, use user provide count_threshold')
}
if (sum(index_count>NB_count_thresh)>1000){
	print('use top 1000 lim')
	NB_count_thresh = index_count[order(-index_count)][1000]
}
print(NB_count_thresh)
if (NB_count_thresh<dim(d1_raw_binary)[2]){
	NB_count_thresh = dim(d1_raw_binary)[2]+1
}
} else if ((dim(d1_raw_binary)[2]+1)<count_threshold){
	NB_count_thresh = count_threshold
} else {
	NB_count_thresh = dim(d1_raw_binary)[2]+1
}



print(NB_count_thresh)
#NB_count_thresh = 50
### filter rare clusters
d1_raw_binary_index_X = rownames(table(d1_raw_binary_index)[table(d1_raw_binary_index)<NB_count_thresh])

### rare cluster to X
for (is in d1_raw_binary_index_X){
	d1_raw_binary_index[d1_raw_binary_index==is] = 'X_X_X_X_X_X'
}
print(length(table(d1_raw_binary_index)))

### get noise
d1_raw_noise = d1_raw+matrix(runif(dim(d1_raw)[1]*dim(d1_raw)[2], min = 0, max = 0.1), ncol=dim(d1_raw)[2])
new_id = d1_raw_binary_index

### QDA clustering
for (i in c(1:500)){
	print(i)
	print('before qda')
	print(table(as.vector(new_id)))
	qda_model = qda(d1_raw_noise, as.character(new_id))
	new_id_new = predict(qda_model, data = d1_raw)$class
	number_of_change = sum(as.character(new_id_new)!=new_id)
	print(number_of_change)
	new_id = as.character(new_id_new)
	print('sum(is.na(new_id))')
	print(sum(is.na(new_id)))
	### get rare
	new_id_index_X = rownames(table(new_id))[table(new_id)<NB_count_thresh]
	for (is in new_id_index_X){
		new_id[new_id==is] = 'X_X_X_X_X_X'
	}
	print('after qda')
	print(new_id_index_X)
	print(table(new_id))
	print(length(table(new_id)))
	if (number_of_change==0){
		break
	}
}

### get mean signal
mean_mat_log2 = c()
mean_mat = c()
number_vec = c()
for (is in rownames(table(new_id))){
	each_cluster_mean_sig = colMeans(d1_raw[new_id==is,])
	mean_mat_log2 = rbind(mean_mat_log2, each_cluster_mean_sig)
	each_cluster_mean_sig = colMeans(d1_raw0[new_id==is,])
	mean_mat = rbind(mean_mat, each_cluster_mean_sig)
	number_vec = c(number_vec, sum(new_id==is))
}

mk_color = read.table('tad.ideasp.colors.tad_and_bd_and_CP.txt', header=F)
mk_color = read.table('tad.ideasp.colors.tadbdonly.txt', header=F)


### get state color
mk_color_vec = t(apply(mk_color, 1, function(x) as.numeric(unlist(strsplit(x, ',')))))
new_state_colors = stateColor(mean_mat, mk_color_vec[-1,])
mk_color_rgb = c()
for (i in 1:dim(mk_color_vec)[1]){
	rgb_i = rgb(mk_color_vec[i,1]/255, mk_color_vec[i,2]/255, mk_color_vec[i,3]/255)
	mk_color_rgb = c(mk_color_rgb, rgb_i)
}
mk_color_rgb = c(mk_color_rgb[-1], "#FFFFFF")

### mean signal color
lim_u = max(mean_mat)
lim_l = min(mean_mat)
color_mat0 = c()
for (i in 1:dim(mean_mat)[1]){
color_row = c()
color_row_code = c()
for (j in 1:dim(mean_mat)[2]){
sig = mean_mat[i,j]
sig_p = (abs(sig-lim_l))/(lim_u-lim_l)
r = ((255-16)*(1-sig_p)+16)/255
g= ((255-78)*(1-sig_p)+78)/255
b = ((255-139)*(1-sig_p)+139)/255
rgb_color = rgb(r,g,b)
color_row[j] = rgb_color
}
color_mat0 = rbind(color_mat0, color_row)
}

color_mat = cbind(color_mat0, new_state_colors[,2])
color_mat_rownames = apply(cbind(rownames(table(new_id)), number_vec, round(number_vec/sum(number_vec), 2)), 1, function(x) paste(x[1], ': ', x[2], ' (', x[3], ')', sep=''))
rownames(color_mat) = c(color_mat_rownames)
colnames(color_mat) = c(2:(dim(color_mat)[2]), 'S')

### hclust
hclust_meansig = hclust(dist(mean_mat))
color_mat_hclust = color_mat[hclust_meansig$order,]
color_mat_hclust = rbind(color_mat_hclust, mk_color_rgb)
rownames(color_mat_hclust) = c(color_mat_rownames,'MK')
color_heatmap(color_mat_hclust, 'qda.cluster.7.statecolor.pdf', pdf, 'gray', 7, 7)


output_newid = cbind(new_id, new_id)
new_id_names = rownames(table(new_id))
for (i in c(1:length(new_id_names))){
color_tmp = new_state_colors[i,1]
output_newid[new_id==new_id_names[i],2] = color_tmp
}

bed = read.table('TAD_regions_tadbdonly.txt', header=F, sep=' ')[,-4]
write.table(cbind(bed, output_newid[,1], rep(1000, dim(bed)[1]), rep('.', dim(bed)[1]), bed[,2], bed[,3], output_newid[,2]), 'output_newid_color.bedgraph', quote=F, col.names=F, row.names=F, sep='\t')


###

#time ~/group/software/ucsc/bedToBigBed output_newid_color.bedgraph ~/group/genome/mm10/mm10.1to19_X.10kb.genome cluster.tad.qda.bb



