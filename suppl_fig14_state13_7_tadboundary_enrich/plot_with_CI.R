
### read state file
HPC7_13_ER4_13 = read.table('HPC7.13.ER4.13.mergedpk.bed', header=F)
HPC7_13_ER4_7 = read.table('HPC7.13.ER4.7.mergedpk.bed', header=F)
HPC7_1_ER4_1 = read.table('HPC7_1.ER4_1.uniq.OnTad.joint.boundary.bed', header=F)
HPC7_1_ER4_0 = read.table('HPC7_1.ER4_0.uniq.OnTad.joint.boundary.bed', header=F)

### get chr list
chr_list = rownames(table(HPC7_13_ER4_13[,1]))
#chr_list = chr_list[1:2]
chr_wg = read.table('~/group/genome/mm10/mm10.1to19_X.genome', header=F)[,2]
exp_win_list = c(0, 10000, 20000, 30000, 40000, 50000, 60000, 70000, 80000, 90000, 100000)
exp_win_list_str = as.character(exp_win_list)
exp_win_list_str[length(exp_win_list_str)] = '100000'


### state chr count
HPC7_13_ER4_13_chrnum = table(HPC7_13_ER4_13[,1])
HPC7_13_ER4_7_chrnum = table(HPC7_13_ER4_7[,1])

### UP
HPC7_1_ER4_0_HPC7_13_ER4_13_enrich_up = c()
HPC7_1_ER4_0_HPC7_13_ER4_7_enrich_up = c()
HPC7_1_ER4_1_HPC7_13_ER4_13_enrich_up = c()
HPC7_1_ER4_1_HPC7_13_ER4_7_enrich_up = c()
for (i in 1:length(exp_win_list)){
expwin_i = exp_win_list_str[i]
HPC7_1_ER4_0_HPC7_13_ER4_13_filenames = paste('HPC7_1.ER4_0.HPC7_13.ER4_13.uniq.OnTad.joint.boundary.up', expwin_i, '.bed', sep='')
HPC7_1_ER4_0_HPC7_13_ER4_13 = read.table(HPC7_1_ER4_0_HPC7_13_ER4_13_filenames, header=F)
HPC7_1_ER4_0_HPC7_13_ER4_7_filenames = paste('HPC7_1.ER4_0.HPC7_13.ER4_7.uniq.OnTad.joint.boundary.up', expwin_i, '.bed', sep='')
HPC7_1_ER4_0_HPC7_13_ER4_7 = read.table(HPC7_1_ER4_0_HPC7_13_ER4_7_filenames, header=F)
HPC7_1_ER4_1_HPC7_13_ER4_13_filenames = paste('HPC7_1.ER4_1.HPC7_13.ER4_13.uniq.OnTad.joint.boundary.up', expwin_i, '.bed', sep='')
HPC7_1_ER4_1_HPC7_13_ER4_13 = read.table(HPC7_1_ER4_1_HPC7_13_ER4_13_filenames, header=F)
HPC7_1_ER4_1_HPC7_13_ER4_7_filenames = paste('HPC7_1.ER4_1.HPC7_13.ER4_7.uniq.OnTad.joint.boundary.up', expwin_i, '.bed', sep='')
HPC7_1_ER4_1_HPC7_13_ER4_7 = read.table(HPC7_1_ER4_1_HPC7_13_ER4_7_filenames, header=F)
### get enrich
print('get enrich')
HPC7_1_ER4_0_HPC7_13_ER4_13_enrich_up_i = c()
HPC7_1_ER4_0_HPC7_13_ER4_7_enrich_up_i = c()
HPC7_1_ER4_1_HPC7_13_ER4_13_enrich_up_i = c()
HPC7_1_ER4_1_HPC7_13_ER4_7_enrich_up_i = c()
for (j in 1:length(chr_list)){
chr_j = chr_list[j]
print(chr_j)
HPC7_13_ER4_13_chrnum_j = HPC7_13_ER4_13_chrnum[j]
HPC7_13_ER4_7_chrnum_j = HPC7_13_ER4_7_chrnum[j]
### exp 
HPC7_1_ER4_0_HPC7_13_ER4_13_exp = (sum(HPC7_1_ER4_0[,1]==chr_j)*10000)/(chr_wg[j]) * (sum(HPC7_13_ER4_13_chrnum_j))
HPC7_1_ER4_0_HPC7_13_ER4_7_exp = (sum(HPC7_1_ER4_0[,1]==chr_j)*10000)/(chr_wg[j]) * (sum(HPC7_13_ER4_7_chrnum_j))
HPC7_1_ER4_1_HPC7_13_ER4_13_exp = (sum(HPC7_1_ER4_1[,1]==chr_j)*10000)/(chr_wg[j]) * (sum(HPC7_13_ER4_13_chrnum_j))
HPC7_1_ER4_1_HPC7_13_ER4_7_exp = (sum(HPC7_1_ER4_1[,1]==chr_j)*10000)/(chr_wg[j]) * (sum(HPC7_13_ER4_7_chrnum_j))
### enrich
HPC7_1_ER4_0_HPC7_13_ER4_13_enrich = sum(HPC7_1_ER4_0_HPC7_13_ER4_13[HPC7_1_ER4_0_HPC7_13_ER4_13[,1]==chr_j,4])/HPC7_1_ER4_0_HPC7_13_ER4_13_exp
HPC7_1_ER4_0_HPC7_13_ER4_7_enrich = sum(HPC7_1_ER4_0_HPC7_13_ER4_7[HPC7_1_ER4_0_HPC7_13_ER4_7[,1]==chr_j,4])/HPC7_1_ER4_0_HPC7_13_ER4_7_exp
HPC7_1_ER4_1_HPC7_13_ER4_13_enrich = sum(HPC7_1_ER4_1_HPC7_13_ER4_13[HPC7_1_ER4_1_HPC7_13_ER4_13[,1]==chr_j,4])/HPC7_1_ER4_1_HPC7_13_ER4_13_exp
HPC7_1_ER4_1_HPC7_13_ER4_7_enrich = sum(HPC7_1_ER4_1_HPC7_13_ER4_7[HPC7_1_ER4_1_HPC7_13_ER4_7[,1]==chr_j,4])/HPC7_1_ER4_1_HPC7_13_ER4_7_exp
### 
HPC7_1_ER4_0_HPC7_13_ER4_13_enrich_up_i = c(HPC7_1_ER4_0_HPC7_13_ER4_13_enrich_up_i, HPC7_1_ER4_0_HPC7_13_ER4_13_enrich)
HPC7_1_ER4_0_HPC7_13_ER4_7_enrich_up_i = c(HPC7_1_ER4_0_HPC7_13_ER4_7_enrich_up_i, HPC7_1_ER4_0_HPC7_13_ER4_7_enrich)
HPC7_1_ER4_1_HPC7_13_ER4_13_enrich_up_i = c(HPC7_1_ER4_1_HPC7_13_ER4_13_enrich_up_i, HPC7_1_ER4_1_HPC7_13_ER4_13_enrich)
HPC7_1_ER4_1_HPC7_13_ER4_7_enrich_up_i = c(HPC7_1_ER4_1_HPC7_13_ER4_7_enrich_up_i, HPC7_1_ER4_1_HPC7_13_ER4_7_enrich)
}
### 
HPC7_1_ER4_0_HPC7_13_ER4_13_enrich_up = cbind(HPC7_1_ER4_0_HPC7_13_ER4_13_enrich_up, HPC7_1_ER4_0_HPC7_13_ER4_13_enrich_up_i)
HPC7_1_ER4_0_HPC7_13_ER4_7_enrich_up = cbind(HPC7_1_ER4_0_HPC7_13_ER4_7_enrich_up, HPC7_1_ER4_0_HPC7_13_ER4_7_enrich_up_i)
HPC7_1_ER4_1_HPC7_13_ER4_13_enrich_up = cbind(HPC7_1_ER4_1_HPC7_13_ER4_13_enrich_up, HPC7_1_ER4_1_HPC7_13_ER4_13_enrich_up_i)
HPC7_1_ER4_1_HPC7_13_ER4_7_enrich_up = cbind(HPC7_1_ER4_1_HPC7_13_ER4_7_enrich_up, HPC7_1_ER4_1_HPC7_13_ER4_7_enrich_up_i)
}
### add names
colnames(HPC7_1_ER4_0_HPC7_13_ER4_13_enrich_up) = exp_win_list
colnames(HPC7_1_ER4_0_HPC7_13_ER4_7_enrich_up) = exp_win_list
colnames(HPC7_1_ER4_1_HPC7_13_ER4_13_enrich_up) = exp_win_list
colnames(HPC7_1_ER4_1_HPC7_13_ER4_7_enrich_up) = exp_win_list
rownames(HPC7_1_ER4_0_HPC7_13_ER4_13_enrich_up) = chr_list
rownames(HPC7_1_ER4_0_HPC7_13_ER4_7_enrich_up) = chr_list
rownames(HPC7_1_ER4_1_HPC7_13_ER4_13_enrich_up) = chr_list
rownames(HPC7_1_ER4_1_HPC7_13_ER4_7_enrich_up) = chr_list


### DOWN
HPC7_1_ER4_0_HPC7_13_ER4_13_enrich_down = c()
HPC7_1_ER4_0_HPC7_13_ER4_7_enrich_down = c()
HPC7_1_ER4_1_HPC7_13_ER4_13_enrich_down = c()
HPC7_1_ER4_1_HPC7_13_ER4_7_enrich_down = c()
for (i in 1:length(exp_win_list)){
expwin_i = exp_win_list_str[i]
HPC7_1_ER4_0_HPC7_13_ER4_13_filenames = paste('HPC7_1.ER4_0.HPC7_13.ER4_13.uniq.OnTad.joint.boundary.down', expwin_i, '.bed', sep='')
HPC7_1_ER4_0_HPC7_13_ER4_13 = read.table(HPC7_1_ER4_0_HPC7_13_ER4_13_filenames, header=F)
HPC7_1_ER4_0_HPC7_13_ER4_7_filenames = paste('HPC7_1.ER4_0.HPC7_13.ER4_7.uniq.OnTad.joint.boundary.down', expwin_i, '.bed', sep='')
HPC7_1_ER4_0_HPC7_13_ER4_7 = read.table(HPC7_1_ER4_0_HPC7_13_ER4_7_filenames, header=F)
HPC7_1_ER4_1_HPC7_13_ER4_13_filenames = paste('HPC7_1.ER4_1.HPC7_13.ER4_13.uniq.OnTad.joint.boundary.down', expwin_i, '.bed', sep='')
HPC7_1_ER4_1_HPC7_13_ER4_13 = read.table(HPC7_1_ER4_1_HPC7_13_ER4_13_filenames, header=F)
HPC7_1_ER4_1_HPC7_13_ER4_7_filenames = paste('HPC7_1.ER4_1.HPC7_13.ER4_7.uniq.OnTad.joint.boundary.down', expwin_i, '.bed', sep='')
HPC7_1_ER4_1_HPC7_13_ER4_7 = read.table(HPC7_1_ER4_1_HPC7_13_ER4_7_filenames, header=F)
### get enrich
print('get enrich')
HPC7_1_ER4_0_HPC7_13_ER4_13_enrich_down_i = c()
HPC7_1_ER4_0_HPC7_13_ER4_7_enrich_down_i = c()
HPC7_1_ER4_1_HPC7_13_ER4_13_enrich_down_i = c()
HPC7_1_ER4_1_HPC7_13_ER4_7_enrich_down_i = c()
for (j in 1:length(chr_list)){
chr_j = chr_list[j]
print(chr_j)
HPC7_13_ER4_13_chrnum_j = HPC7_13_ER4_13_chrnum[j]
HPC7_13_ER4_7_chrnum_j = HPC7_13_ER4_7_chrnum[j]
### exp 
HPC7_1_ER4_0_HPC7_13_ER4_13_exp = (sum(HPC7_1_ER4_0[,1]==chr_j)*10000)/(chr_wg[j]) * (sum(HPC7_13_ER4_13_chrnum_j))
HPC7_1_ER4_0_HPC7_13_ER4_7_exp = (sum(HPC7_1_ER4_0[,1]==chr_j)*10000)/(chr_wg[j]) * (sum(HPC7_13_ER4_7_chrnum_j))
HPC7_1_ER4_1_HPC7_13_ER4_13_exp = (sum(HPC7_1_ER4_1[,1]==chr_j)*10000)/(chr_wg[j]) * (sum(HPC7_13_ER4_13_chrnum_j))
HPC7_1_ER4_1_HPC7_13_ER4_7_exp = (sum(HPC7_1_ER4_1[,1]==chr_j)*10000)/(chr_wg[j]) * (sum(HPC7_13_ER4_7_chrnum_j))
### enrich
HPC7_1_ER4_0_HPC7_13_ER4_13_enrich = sum(HPC7_1_ER4_0_HPC7_13_ER4_13[HPC7_1_ER4_0_HPC7_13_ER4_13[,1]==chr_j,4])/HPC7_1_ER4_0_HPC7_13_ER4_13_exp
HPC7_1_ER4_0_HPC7_13_ER4_7_enrich = sum(HPC7_1_ER4_0_HPC7_13_ER4_7[HPC7_1_ER4_0_HPC7_13_ER4_7[,1]==chr_j,4])/HPC7_1_ER4_0_HPC7_13_ER4_7_exp
HPC7_1_ER4_1_HPC7_13_ER4_13_enrich = sum(HPC7_1_ER4_1_HPC7_13_ER4_13[HPC7_1_ER4_1_HPC7_13_ER4_13[,1]==chr_j,4])/HPC7_1_ER4_1_HPC7_13_ER4_13_exp
HPC7_1_ER4_1_HPC7_13_ER4_7_enrich = sum(HPC7_1_ER4_1_HPC7_13_ER4_7[HPC7_1_ER4_1_HPC7_13_ER4_7[,1]==chr_j,4])/HPC7_1_ER4_1_HPC7_13_ER4_7_exp
### 
HPC7_1_ER4_0_HPC7_13_ER4_13_enrich_down_i = c(HPC7_1_ER4_0_HPC7_13_ER4_13_enrich_down_i, HPC7_1_ER4_0_HPC7_13_ER4_13_enrich)
HPC7_1_ER4_0_HPC7_13_ER4_7_enrich_down_i = c(HPC7_1_ER4_0_HPC7_13_ER4_7_enrich_down_i, HPC7_1_ER4_0_HPC7_13_ER4_7_enrich)
HPC7_1_ER4_1_HPC7_13_ER4_13_enrich_down_i = c(HPC7_1_ER4_1_HPC7_13_ER4_13_enrich_down_i, HPC7_1_ER4_1_HPC7_13_ER4_13_enrich)
HPC7_1_ER4_1_HPC7_13_ER4_7_enrich_down_i = c(HPC7_1_ER4_1_HPC7_13_ER4_7_enrich_down_i, HPC7_1_ER4_1_HPC7_13_ER4_7_enrich)
}
### 
HPC7_1_ER4_0_HPC7_13_ER4_13_enrich_down = cbind(HPC7_1_ER4_0_HPC7_13_ER4_13_enrich_down, HPC7_1_ER4_0_HPC7_13_ER4_13_enrich_down_i)
HPC7_1_ER4_0_HPC7_13_ER4_7_enrich_down = cbind(HPC7_1_ER4_0_HPC7_13_ER4_7_enrich_down, HPC7_1_ER4_0_HPC7_13_ER4_7_enrich_down_i)
HPC7_1_ER4_1_HPC7_13_ER4_13_enrich_down = cbind(HPC7_1_ER4_1_HPC7_13_ER4_13_enrich_down, HPC7_1_ER4_1_HPC7_13_ER4_13_enrich_down_i)
HPC7_1_ER4_1_HPC7_13_ER4_7_enrich_down = cbind(HPC7_1_ER4_1_HPC7_13_ER4_7_enrich_down, HPC7_1_ER4_1_HPC7_13_ER4_7_enrich_down_i)
}
### add names
colnames(HPC7_1_ER4_0_HPC7_13_ER4_13_enrich_down) = exp_win_list
colnames(HPC7_1_ER4_0_HPC7_13_ER4_7_enrich_down) = exp_win_list
colnames(HPC7_1_ER4_1_HPC7_13_ER4_13_enrich_down) = exp_win_list
colnames(HPC7_1_ER4_1_HPC7_13_ER4_7_enrich_down) = exp_win_list
rownames(HPC7_1_ER4_0_HPC7_13_ER4_13_enrich_down) = chr_list
rownames(HPC7_1_ER4_0_HPC7_13_ER4_7_enrich_down) = chr_list
rownames(HPC7_1_ER4_1_HPC7_13_ER4_13_enrich_down) = chr_list
rownames(HPC7_1_ER4_1_HPC7_13_ER4_7_enrich_down) = chr_list


### matrix
HPC7_1_ER4_0_HPC7_13_ER4_13_enrich_updown = cbind(HPC7_1_ER4_0_HPC7_13_ER4_13_enrich_up[,dim(HPC7_1_ER4_0_HPC7_13_ER4_13_enrich_up)[2]:2], HPC7_1_ER4_0_HPC7_13_ER4_13_enrich_down)
HPC7_1_ER4_0_HPC7_13_ER4_7_enrich_updown = cbind(HPC7_1_ER4_0_HPC7_13_ER4_7_enrich_up[,dim(HPC7_1_ER4_0_HPC7_13_ER4_13_enrich_up)[2]:2], HPC7_1_ER4_0_HPC7_13_ER4_7_enrich_down)
HPC7_1_ER4_1_HPC7_13_ER4_13_enrich_updown = cbind(HPC7_1_ER4_1_HPC7_13_ER4_13_enrich_up[,dim(HPC7_1_ER4_0_HPC7_13_ER4_13_enrich_up)[2]:2], HPC7_1_ER4_1_HPC7_13_ER4_13_enrich_down)
HPC7_1_ER4_1_HPC7_13_ER4_7_enrich_updown = cbind(HPC7_1_ER4_1_HPC7_13_ER4_7_enrich_up[,dim(HPC7_1_ER4_0_HPC7_13_ER4_13_enrich_up)[2]:2], HPC7_1_ER4_1_HPC7_13_ER4_7_enrich_down)

### ggplot2
require(ggplot2)


ci = function(x){
x_ci1 = sd(x)/sqrt(length(x))*1.96
return(x_ci1)
}

### get df
exp_win_list_updown = c(-exp_win_list[length(exp_win_list):2], exp_win_list)
i=1
HPC7_1_ER4_0_enrich_dfa = as.data.frame(cbind(as.data.frame(mean(as.numeric(HPC7_1_ER4_0_HPC7_13_ER4_13_enrich_updown[,i]))), exp_win_list_updown[i]/1000, 'HPC7_13_ER4_13', as.numeric(ci(HPC7_1_ER4_0_HPC7_13_ER4_13_enrich_updown[,i])) ))
for (i in 2:dim(HPC7_1_ER4_0_HPC7_13_ER4_13_enrich_updown)[2]){
HPC7_1_ER4_0_enrich_df_i = as.data.frame(cbind(as.data.frame(mean(as.numeric(HPC7_1_ER4_0_HPC7_13_ER4_13_enrich_updown[,i]))), exp_win_list_updown[i]/1000, 'HPC7_13_ER4_13', as.numeric(ci(HPC7_1_ER4_0_HPC7_13_ER4_13_enrich_updown[,i])) ))
HPC7_1_ER4_0_enrich_dfa = rbind(HPC7_1_ER4_0_enrich_dfa, HPC7_1_ER4_0_enrich_df_i)
}

i=1
HPC7_1_ER4_0_enrich_dfb = as.data.frame(cbind(as.data.frame(mean(as.numeric(HPC7_1_ER4_0_HPC7_13_ER4_7_enrich_updown[,i]))), exp_win_list_updown[i]/1000, 'HPC7_13_ER4_7', as.numeric(ci(HPC7_1_ER4_0_HPC7_13_ER4_7_enrich_updown[,i])) ))
for (i in 2:dim(HPC7_1_ER4_0_HPC7_13_ER4_7_enrich_updown)[2]){
HPC7_1_ER4_0_enrich_df_i = as.data.frame(cbind(as.data.frame(mean(as.numeric(HPC7_1_ER4_0_HPC7_13_ER4_7_enrich_updown[,i]))), exp_win_list_updown[i]/1000, 'HPC7_13_ER4_7', as.numeric(ci(HPC7_1_ER4_0_HPC7_13_ER4_7_enrich_updown[,i])) ))
HPC7_1_ER4_0_enrich_dfb = rbind(HPC7_1_ER4_0_enrich_dfb, HPC7_1_ER4_0_enrich_df_i)
}
colnames(HPC7_1_ER4_0_enrich_dfa) = c('Enrich', 'win', 'TADboundary', 'ci')
colnames(HPC7_1_ER4_0_enrich_dfb) = c('Enrich', 'win', 'TADboundary', 'ci')
HPC7_1_ER4_0_enrich_df = rbind(HPC7_1_ER4_0_enrich_dfa, HPC7_1_ER4_0_enrich_dfb)

pd <- position_dodge(0.1)
pdf('HPC7_1_ER4_0_enrich.pdf', width=10, height=4)
ggplot(HPC7_1_ER4_0_enrich_df, aes(x=win, y=Enrich, colour=TADboundary)) + 
    geom_errorbar(aes(ymin=Enrich-ci, ymax=Enrich+ci), width=2, position=pd) +
    geom_line(position=pd) +
    geom_point(position=pd) + 
    theme(text = element_text(size=20))
dev.off()

### get df
exp_win_list_updown = c(-exp_win_list[length(exp_win_list):2], exp_win_list)
i=1
HPC7_1_ER4_1_enrich_dfa = as.data.frame(cbind(as.data.frame(mean(as.numeric(HPC7_1_ER4_1_HPC7_13_ER4_13_enrich_updown[,i]))), exp_win_list_updown[i]/1000, 'HPC7_13_ER4_13', as.numeric(ci(HPC7_1_ER4_1_HPC7_13_ER4_13_enrich_updown[,i])) ))
for (i in 2:dim(HPC7_1_ER4_1_HPC7_13_ER4_13_enrich_updown)[2]){
HPC7_1_ER4_1_enrich_df_i = as.data.frame(cbind(as.data.frame(mean(as.numeric(HPC7_1_ER4_1_HPC7_13_ER4_13_enrich_updown[,i]))), exp_win_list_updown[i]/1000, 'HPC7_13_ER4_13', as.numeric(ci(HPC7_1_ER4_1_HPC7_13_ER4_13_enrich_updown[,i])) ))
HPC7_1_ER4_1_enrich_dfa = rbind(HPC7_1_ER4_1_enrich_dfa, HPC7_1_ER4_1_enrich_df_i)
}

i=1
HPC7_1_ER4_1_enrich_dfb = as.data.frame(cbind(as.data.frame(mean(as.numeric(HPC7_1_ER4_1_HPC7_13_ER4_7_enrich_updown[,i]))), exp_win_list_updown[i]/1000, 'HPC7_13_ER4_7', as.numeric(ci(HPC7_1_ER4_1_HPC7_13_ER4_7_enrich_updown[,i])) ))
for (i in 2:dim(HPC7_1_ER4_1_HPC7_13_ER4_7_enrich_updown)[2]){
HPC7_1_ER4_1_enrich_df_i = as.data.frame(cbind(as.data.frame(mean(as.numeric(HPC7_1_ER4_1_HPC7_13_ER4_7_enrich_updown[,i]))), exp_win_list_updown[i]/1000, 'HPC7_13_ER4_7', as.numeric(ci(HPC7_1_ER4_1_HPC7_13_ER4_7_enrich_updown[,i])) ))
HPC7_1_ER4_1_enrich_dfb = rbind(HPC7_1_ER4_1_enrich_dfb, HPC7_1_ER4_1_enrich_df_i)
}
colnames(HPC7_1_ER4_1_enrich_dfa) = c('Enrich', 'win', 'TADboundary', 'ci')
colnames(HPC7_1_ER4_1_enrich_dfb) = c('Enrich', 'win', 'TADboundary', 'ci')
HPC7_1_ER4_1_enrich_df = rbind(HPC7_1_ER4_1_enrich_dfa, HPC7_1_ER4_1_enrich_dfb)

pd <- position_dodge(0.1)
pdf('HPC7_1_ER4_1_enrich.pdf', width=10, height=4)
ggplot(HPC7_1_ER4_1_enrich_df, aes(x=win, y=Enrich, colour=TADboundary)) + 
    geom_errorbar(aes(ymin=Enrich-ci, ymax=Enrich+ci), width=2, position=pd) +
    geom_line(position=pd) +
    geom_point(position=pd) +
    theme(text = element_text(size=20))
dev.off()





