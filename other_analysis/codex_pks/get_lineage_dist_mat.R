library(pheatmap)

IS_vec = read.table('IS_mat.txt', header=F)
sig = read.table('IS_sig_mat.txt', header=F)
ct_list = read.table('IS_index_ct.txt', header=F)

### distance seperate
get_c1_c2_vardist = function(sig_i, c1, c2){
### get start cell pairs
sig_d1 = c()
for (ct in c1){
sig_d1_tmp = sig_i[,ct_list[ct_list[,2]==ct,1]]
sig_d1 = cbind(sig_d1, sig_d1_tmp)
}
sig_d2 = c()
for (ct in c2){
sig_d2_tmp = sig_i[,ct_list[ct_list[,2]==ct,1]]
sig_d2 = cbind(sig_d2, sig_d2_tmp)
}
sig_d12 = cbind(sig_d1, sig_d2)
### get variance within branch
#var_1 = sum((sig_d1-rowMeans(sig_d1))^2)/length(c1)
#var_2 = sum((sig_d2-rowMeans(sig_d2))^2)/length(c2)
#var_12 = sum((sig_d12-rowMeans(sig_d12))^2)/(length(c1)+length(c2))
var_1 = sum(apply(sig_d1, 1, function(x) sd(x)/mean(x)))
var_2 = sum(apply(sig_d2, 1, function(x) sd(x)/mean(x)))
var_12 = sum(apply(sig_d12, 1, function(x) sd(x)/mean(x)))
### get 1-var1_ratio
small_num = 0
var1_ratio = (var_1+small_num)/(var_12+small_num)
var2_ratio = (var_2+small_num)/(var_12+small_num)
var12_ratio = (var_1+var_2+small_num)/(var_12+small_num)
return(c(var1_ratio, var2_ratio, var12_ratio))
}

### distance seperate
get_c1_c2_vardist = function(sig_i, c1, c2){
### get start cell pairs
sig_d1 = c()
for (ct in c1){
sig_d1_tmp = sig_i[,ct_list[ct_list[,2]==ct,1]]
sig_d1 = cbind(sig_d1, sig_d1_tmp)
}
if (dim(sig_d1)[2]>2){
sig_d1 = cbind(sig_d1[,1], rowMeans(sig_d1[,-1]))
}
sig_d2 = c()
for (ct in c2){
sig_d2_tmp = sig_i[,ct_list[ct_list[,2]==ct,1]]
sig_d2 = cbind(sig_d2, sig_d2_tmp)
}
if (dim(sig_d1)[2]>2){
sig_d2 = cbind(sig_d2[,1], rowMeans(sig_d2[,-1]))
}
sig_d12 = cbind(sig_d1, sig_d2)
### get variance within branch
#var_1 = sum((sig_d1-rowMeans(sig_d1))^2)/length(c1)
#var_2 = sum((sig_d2-rowMeans(sig_d2))^2)/length(c2)
#var_12 = sum((sig_d12-rowMeans(sig_d12))^2)/(length(c1)+length(c2))
var_1 = sum(apply(sig_d1, 1, function(x) sd(x)/mean(x)))
var_2 = sum(apply(sig_d2, 1, function(x) sd(x)/mean(x)))
var_12 = sum(apply(sig_d12, 1, function(x) sd(x)/mean(x)))
### get 1-var1_ratio
small_num = 0
var1_ratio = (var_1+small_num)/(var_12+small_num)
var2_ratio = (var_2+small_num)/(var_12+small_num)
var12_ratio = (var_1+var_2+small_num)/(var_12+small_num)
return(c(var1_ratio, var2_ratio, var12_ratio))
}

get_all_IS_vardist = function(sig, IS_vec, c1_1, c2_1){
p1_mat = c()
for (i in c(1:dim(IS_vec)[1])){
IS_i = IS_vec[i,1]
#IS_i = IS_vec[1,1]
### get IS_i signal matrix
sig_i = sig[sig[,2]==IS_i,-c(1,2)]
p1 = get_c1_c2_vardist(sig_i, c1_1, c2_1)
p1_mat = rbind(p1_mat, p1)
}
return(p1_mat)
}

c1_1 = c('CFUE')
c2_1 = c('CFUMK')
p1_mat = get_all_IS_vardist(sig, IS_vec, c1_1, c2_1)

c1_1 = c('CFUE', 'ERY', 'ERY_fl')
c2_1 = c('CFUMK', 'iMK')
p2_mat = get_all_IS_vardist(sig, IS_vec, c1_1, c2_1)

c1_1 = c('MEP')
c2_1 = c('GMP')
p3_mat = get_all_IS_vardist(sig, IS_vec, c1_1, c2_1)

c1_1 = c('MEP', 'CFUE', 'ERY', 'ERY_fl', 'CFUMK', 'iMK')
c2_1 = c('GMP', 'NEU', 'MON')
p4_mat = get_all_IS_vardist(sig, IS_vec, c1_1, c2_1)

c1_1 = c('CMP', 'MEP', 'CFUE', 'ERY', 'ERY_fl', 'CFUMK', 'iMK', 'GMP', 'NEU', 'MON')
c2_1 = c('B', 'NK', 'T_CD4', 'T_CD8')
p5_mat = get_all_IS_vardist(sig, IS_vec, c1_1, c2_1)


p_all = cbind(p2_mat[,], p4_mat[,], p5_mat[,])

corlo_lim_scale_max = max(p_all)
corlo_lim_scale_min = min(p_all)
breaksList_scale = seq(corlo_lim_scale_min, corlo_lim_scale_max, by = 0.001)
my_colorbar_scale=colorRampPalette(c('red', 'white'))(n = length(breaksList_scale))
#JI_dist = dist(t(log10(JI_mat+1e-4)))
pdf('branch_vardist.pheatmap.pdf', width = 30, height = 30)
pheatmap(as.matrix(p_all), color=my_colorbar_scale, breaks = breaksList_scale, annotation_names_col = FALSE, cluster_cols=FALSE, cluster_rows=FALSE)
dev.off()

p_all = cbind(p2_mat[,3], p4_mat[,3], p5_mat[,3])

corlo_lim_scale_max = max(p_all)
corlo_lim_scale_min = min(p_all)
breaksList_scale = seq(corlo_lim_scale_min, corlo_lim_scale_max, by = 0.001)
my_colorbar_scale=colorRampPalette(c('red', 'white'))(n = length(breaksList_scale))
#JI_dist = dist(t(log10(JI_mat+1e-4)))
pdf('branch_vardist.pheatmap.3.pdf', width = 30, height = 30)
pheatmap(as.matrix(p_all), color=my_colorbar_scale, breaks = breaksList_scale, annotation_names_col = FALSE, cluster_cols=FALSE, cluster_rows=FALSE)
dev.off()



get_all_IS_vardist = function(sig, IS_vec, c1_1, c2_1){
p1_mat = c()
for (i in c(1:dim(IS_vec)[1])){
IS_i = IS_vec[i,1]
#IS_i = IS_vec[1,1]
### get IS_i signal matrix
sig_i = sig[sig[,2]==IS_i,1:2]
sig_i = t(apply(sig_i, 1, function(x) as.numeric(unlist(strsplit(x[2], '_')))))
#print(head(sig_i))
p1 = get_c1_c2_vardist(sig_i, c1_1, c2_1)
p1_mat = rbind(p1_mat, p1)
}
return(p1_mat)
}

### distance seperate
get_c1_c2_vardist = function(sig_i, c1, c2){
### get start cell pairs
sig_d1 = c()
for (ct in c1){
sig_d1_tmp = sig_i[,ct_list[ct_list[,2]==ct,1]]
sig_d1 = cbind(sig_d1, sig_d1_tmp)
}
sig_d2 = c()
for (ct in c2){
sig_d2_tmp = sig_i[,ct_list[ct_list[,2]==ct,1]]
sig_d2 = cbind(sig_d2, sig_d2_tmp)
}
sig_d12 = cbind(sig_d1, sig_d2)
### get variance within branch
#var_1 = sum((sig_d1-rowMeans(sig_d1))^2)/length(c1)
#var_2 = sum((sig_d2-rowMeans(sig_d2))^2)/length(c2)
#var_12 = sum((sig_d12-rowMeans(sig_d12))^2)/(length(c1)+length(c2))
small_num = 1
var_1 = sum(apply(sig_d1, 1, function(x) (sd(x)+small_num)/(mean(x)+small_num)))
var_2 = sum(apply(sig_d2, 1, function(x) (sd(x)+small_num)/(mean(x)+small_num)))
var_12 = sum(apply(sig_d12, 1, function(x) (sd(x)+small_num)/(mean(x)+small_num)))
### get 1-var1_ratio
small_num = 0
var1_ratio = (var_1+small_num)/(var_12+small_num)
var2_ratio = (var_2+small_num)/(var_12+small_num)
var12_ratio = (var_1+var_2+small_num)/(var_12+small_num)
return(c(var1_ratio, var2_ratio, var12_ratio))
}


c1_1 = c('CFUE', 'ERY', 'ERY_fl')
c2_1 = c('CFUMK', 'iMK')
p2_mat = get_all_IS_vardist(sig, IS_vec, c1_1, c2_1)

c1_1 = c('MEP', 'CFUE', 'ERY', 'ERY_fl', 'CFUMK', 'iMK')
c2_1 = c('GMP', 'NEU', 'MON')
p4_mat = get_all_IS_vardist(sig, IS_vec, c1_1, c2_1)

c1_1 = c('CMP', 'MEP', 'CFUE', 'ERY', 'ERY_fl', 'CFUMK', 'iMK', 'GMP', 'NEU', 'MON')
c2_1 = c('B', 'NK', 'T_CD4', 'T_CD8')
p5_mat = get_all_IS_vardist(sig, IS_vec, c1_1, c2_1)


p_all = cbind(p2_mat[,], p4_mat[,], p5_mat[,])

corlo_lim_scale_max = max(p_all[!is.na(p_all)])
corlo_lim_scale_min = min(p_all[!is.na(p_all)])
breaksList_scale = seq(corlo_lim_scale_min, corlo_lim_scale_max, by = 0.001)
my_colorbar_scale=colorRampPalette(c('red', 'white'))(n = length(breaksList_scale))
#JI_dist = dist(t(log10(JI_mat+1e-4)))
pdf('branch_vardist.pheatmap.binary.pdf', width = 30, height = 30)
pheatmap(as.matrix(p_all), color=my_colorbar_scale, breaks = breaksList_scale, annotation_names_col = FALSE, cluster_cols=FALSE, cluster_rows=FALSE)
dev.off()

p_all = cbind(p2_mat[,3], p4_mat[,3], p5_mat[,3])

corlo_lim_scale_max = max(p_all[!is.na(p_all)])
corlo_lim_scale_min = min(p_all[!is.na(p_all)])
breaksList_scale = seq(corlo_lim_scale_min, corlo_lim_scale_max, by = 0.001)
my_colorbar_scale=colorRampPalette(c('red', 'white'))(n = length(breaksList_scale))
#JI_dist = dist(t(log10(JI_mat+1e-4)))
pdf('branch_vardist.pheatmap.3.binary.pdf', width = 30, height = 30)
pheatmap(as.matrix(p_all), color=my_colorbar_scale, breaks = breaksList_scale, annotation_names_col = FALSE, cluster_cols=FALSE, cluster_rows=FALSE)
dev.off()


### distance seperate
get_c1_c2_vardist = function(sig_i, c1, c2){
### get start cell pairs
sig_d1 = c()
for (ct in c1){
sig_d1_tmp = sig_i[,ct_list[ct_list[,2]==ct,1]]
sig_d1 = cbind(sig_d1, sig_d1_tmp)
}
if (dim(sig_d1)[2]>2){
sig_d1 = cbind(sig_d1[,1], rowMeans(sig_d1[,-1]))
}
sig_d2 = c()
for (ct in c2){
sig_d2_tmp = sig_i[,ct_list[ct_list[,2]==ct,1]]
sig_d2 = cbind(sig_d2, sig_d2_tmp)
}
if (dim(sig_d2)[2]>2){
sig_d2 = cbind(sig_d2[,1], rowMeans(sig_d2[,-1]))
}
sig_d12 = cbind(sig_d1, sig_d2)
### get variance within branch
#var_1 = sum((sig_d1-rowMeans(sig_d1))^2)/length(c1)
#var_2 = sum((sig_d2-rowMeans(sig_d2))^2)/length(c2)
#var_12 = sum((sig_d12-rowMeans(sig_d12))^2)/(length(c1)+length(c2))
small_num = 1
var_1 = sum(apply(sig_d1, 1, function(x) (sd(x)+small_num)/(mean(x)+small_num)))
var_2 = sum(apply(sig_d2, 1, function(x) (sd(x)+small_num)/(mean(x)+small_num)))
var_12 = sum(apply(sig_d12, 1, function(x) (sd(x)+small_num)/(mean(x)+small_num)))
### get 1-var1_ratio
small_num = 0
var1_ratio = (var_1+small_num)/(var_12+small_num)
var2_ratio = (var_2+small_num)/(var_12+small_num)
var12_ratio = (var_1+var_2+small_num)/(var_12+small_num)
return(c(var1_ratio, var2_ratio, var12_ratio))
}

get_all_IS_vardist = function(sig, IS_vec, c1_1, c2_1){
p1_mat = c()
for (i in c(1:dim(IS_vec)[1])){
IS_i = IS_vec[i,1]
#IS_i = IS_vec[1,1]
### get IS_i signal matrix
sig_i = sig[sig[,2]==IS_i,1:2]
sig_i = t(apply(sig_i, 1, function(x) as.numeric(unlist(strsplit(x[2], '_')))))
#print(head(sig_i))
p1 = get_c1_c2_vardist(sig_i, c1_1, c2_1)
p1_mat = rbind(p1_mat, p1)
}
return(p1_mat)
}

c1_1 = c('CFUE', 'ERY', 'ERY_fl')
c2_1 = c('CFUMK', 'iMK')
p2_mat = get_all_IS_vardist(sig, IS_vec, c1_1, c2_1)

c1_1 = c('MEP', 'CFUE', 'ERY', 'ERY_fl', 'CFUMK', 'iMK')
c2_1 = c('GMP', 'NEU', 'MON')
p4_mat = get_all_IS_vardist(sig, IS_vec, c1_1, c2_1)

c1_1 = c('CMP', 'MEP', 'CFUE', 'ERY', 'ERY_fl', 'CFUMK', 'iMK', 'GMP', 'NEU', 'MON')
c2_1 = c('B', 'NK', 'T_CD4', 'T_CD8')
p5_mat = get_all_IS_vardist(sig, IS_vec, c1_1, c2_1)


p_all = cbind(p2_mat[,], p4_mat[,], p5_mat[,])

corlo_lim_scale_max = max(p_all[!is.na(p_all)])
corlo_lim_scale_min = min(p_all[!is.na(p_all)])
breaksList_scale = seq(corlo_lim_scale_min, corlo_lim_scale_max, by = 0.001)
my_colorbar_scale=colorRampPalette(c('red', 'white'))(n = length(breaksList_scale))
#JI_dist = dist(t(log10(JI_mat+1e-4)))
pdf('branch_vardist.pheatmap.binary.2.pdf', width = 30, height = 30)
pheatmap(as.matrix(p_all), color=my_colorbar_scale, breaks = breaksList_scale, annotation_names_col = FALSE, cluster_cols=FALSE, cluster_rows=FALSE)
dev.off()

p_all = cbind(p2_mat[,3], p4_mat[,3], p5_mat[,3])

corlo_lim_scale_max = max(p_all[!is.na(p_all)])
corlo_lim_scale_min = min(p_all[!is.na(p_all)])
breaksList_scale = seq(corlo_lim_scale_min, corlo_lim_scale_max, by = 0.001)
my_colorbar_scale=colorRampPalette(c('red', 'white'))(n = length(breaksList_scale))
#JI_dist = dist(t(log10(JI_mat+1e-4)))
pdf('branch_vardist.pheatmap.3.binary.2.pdf', width = 30, height = 30)
pheatmap(as.matrix(p_all), color=my_colorbar_scale, breaks = breaksList_scale, annotation_names_col = FALSE, cluster_cols=FALSE, cluster_rows=FALSE)
dev.off()


   V1     V2
1   1    LSK
2   2   HPC7
3   3    CMP
4   4    MEP
5   5    G1E
6   6    ER4
7   7   CFUE
8   8    ERY
9   9 ERY_fl
10 10  CFUMK
11 11    iMK
12 12    GMP
13 13    MON
14 14    NEU
15 15     NK
16 16      B
17 17  T_CD4
18 18  T_CD8