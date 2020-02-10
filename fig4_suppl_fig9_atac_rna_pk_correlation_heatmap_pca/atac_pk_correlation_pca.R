set.seed(2018)
#signal_matrix = read.table('snapshot20_reproduce_2_16lim/atac_20cell.signal.matrix.no0.txt', header=F)[,c(-1,-2,-3,-4)]
signal_matrix = read.table('atac_20cell.signal.matrix.no0.txt', header=F)[,c(-1,-2,-3,-4)]

names = read.table('signal_list.txt', header=F)[,2]

colnames(signal_matrix) = names

cor_method = 'spearman'
ct_dist = dist(t(signal_matrix))
ct_cor = cor(signal_matrix, method=cor_method)

hc = hclust(ct_dist, method = "average")

library(pheatmap)
my_colorbar=colorRampPalette(c('white', 'red'))(n = 128)
col_breaks = c(seq(0, 2000,length=33))

cor_matrix = cor(signal_matrix, method = cor_method)
cr=round(cor_matrix, digits=2)
library("gplots")
library("RColorBrewer")
my_palette <- c("#1000FFFF","#0047FFFF","#00C8FFFF","#00FFE0FF","#00FF89FF","#00FF33FF","#99FF00FF","#AFFF00FF","#FFD800FF","#FFB800FF","#FFAD00FF","#FF8100FF","#FF5600FF","#FF2B00FF","#FF0000FF")
breaks = c(seq(-.1,1,length=16))

dist.pear <- function(x) as.dist(1-x)
pdf('hc_heatmap.pdf')
#pheatmap(cor_matrix, color=my_colorbar, cluster_cols = TRUE,cluster_rows=TRUE,annotation_names_row=TRUE,annotation_names_col=TRUE,show_rownames=TRUE,show_colnames=TRUE,clustering_method = 'complete')
heatmap.2(cor_matrix, dendrogram="row", col=my_palette, breaks=breaks, trace="none", notecol="black", revC=TRUE, notecex=.5, symm=TRUE, distfun=dist.pear, symkey=FALSE, margins = c(8,8))
dev.off()

pca=prcomp(t(signal_matrix), center = FALSE, scale. = FALSE)
scores=data.frame(pca$x)
pk_loading = pca$rotation

scores_pc1 = as.matrix(-scores[,1])
rownames(scores_pc1) = rownames(scores)
colnames(scores_pc1) = 'PC1'
write.table(scores_pc1, 'negative_pc1_score.txt', quote=FALSE, sep='\t', col.names=FALSE)
scores_pc1_rep1 = as.matrix(-scores[c(2, 4, 6, 8, 10, 13, 15, 17, 19, 21, 23, 25, 27),1])
rownames(scores_pc1_rep1) = rownames(scores)[c(2, 4, 6, 8, 10, 13, 15, 17, 19, 21, 23, 25, 27)]
colnames(scores_pc1_rep1) = 'PC1'
write.table(scores_pc1_rep1, 'negative_pc1_score.rep1.txt', quote=FALSE, sep='\t', col.names=FALSE)
scores_pc1_rep2 = as.matrix(-scores[c(3, 5, 7, 9, 11, 14, 16, 18, 20, 22, 24, 26, 28),1])
rownames(scores_pc1_rep2) = rownames(scores)[c(3, 5, 7, 9, 11, 14, 16, 18, 20, 22, 24, 26, 28)]
colnames(scores_pc1_rep2) = 'PC1'
write.table(scores_pc1_rep2, 'negative_pc1_score.rep2.txt', quote=FALSE, sep='\t', col.names=FALSE)
scores_pc1_0rep = as.matrix(-scores[c(1,12,29,30,31),1])
rownames(scores_pc1_0rep) = rownames(scores)[c(1,12,29,30,31)]
colnames(scores_pc1_0rep) = 'PC1'
write.table(scores_pc1_0rep, 'negative_pc1_score.0rep.txt', quote=FALSE, sep='\t', col.names=FALSE)



contribution_pc1 = apply(signal_matrix, 1, function(x) cor(x, scores$PC1))
contribution_pc2 = apply(signal_matrix, 1, function(x) cor(x, scores$PC2))
contribution_pc3 = apply(signal_matrix, 1, function(x) cor(x, scores$PC3))

pdf('pca_cor.pdf', width=15, height=5)
par(mfrow=c(1,3))
hist(contribution_pc1, breaks=100, xlim=c(-1,1))
box()
hist(contribution_pc2, breaks=100, xlim=c(-1,1))
box()
hist(contribution_pc3, breaks=100, xlim=c(-1,1))
box()
dev.off()

library(pheatmap)
pk_loading_sample = pk_loading[sample(dim(pk_loading)[1], 10000),c(1,2,3)]
pdf('pca_loading.pdf', width=10, height=15)
pheatmap(pk_loading_sample,cluster_cols = FALSE)
dev.off()



cols=ncol(scores)
cells=rownames(scores)#gsub("_.*", "", rownames(scores))
scores=data.frame(scores,cells)
colnames(scores)[cols+1]="cell"

pdf('hc_tree.pdf')
plot(hclust(as.dist(1-cor_matrix), method = "ward.D", members = NULL))
dev.off()



###### PCA analysis
library(ggplot2)
### for cell type
my_palette = c("#8B1C62","#ff3030","#8b7355","#ee7600","#cd5555","#FF0000","#ee2c2c","#ee6363","#008b8b","#ff7f00","#228b22","#cd6600","#8b5a2b","#1874cd","#4f94cd","#8b008b","#7a378b","#68228b")
### for sample
my_palette = c('#8B1C62', '#ff3030', '#ff3030', '#8b7355', '#8b7355', '#ee7600', '#ee7600', '#cd5555', '#cd5555', '#FF0000', '#FF0000', '#ee2c2c', '#ee6363', '#ee6363', '#008b8b', '#008b8b', '#ff7f00', '#ff7f00', '#228b22', '#228b22', '#cd6600', '#cd6600', '#8b5a2b', '#8b5a2b', '#1874cd', '#1874cd', '#4f94cd', '#4f94cd', '#8b008b', '#7a378b', '#68228b')


pdf('pca_pc1_pc2.pdf', useDingbats=FALSE)
cell = colnames(signal_matrix)
ggplot(scores, aes(x=PC2, y=PC1))+geom_point(colour=my_palette, size=4, pch=20)+geom_text(aes(label=rownames(scores)), size=3)+scale_colour_manual(values = my_palette) + theme(legend.position="none")
dev.off()

pdf('pca_pc1_pc3.pdf', useDingbats=FALSE)
cell = colnames(signal_matrix)
ggplot(scores, aes(x=PC3, y=PC1))+geom_point(colour=my_palette, size=4, pch=20)+geom_text(aes(label=rownames(scores)), size=3)+scale_colour_manual(values = my_palette) + theme(legend.position="none")
dev.off()


pdf('pca_pc2_pc3.pdf', useDingbats=FALSE)
cell = colnames(signal_matrix)
ggplot(scores, aes(x=PC2, y=PC3))+geom_point(colour=my_palette, size=4, pch=20)+geom_text(aes(label=rownames(scores)), size=2)+scale_colour_manual(values = my_palette) + theme(legend.position="none")
dev.off()


pdf('pca_pc1.pdf', useDingbats=FALSE, width = 3, height = 7)
cell = colnames(signal_matrix)
ggplot(scores, aes(x=0, y=PC1), col=my_palette)+geom_point(aes(color=cell), size=4, pch=20)+geom_text(aes(label=rownames(scores)), size=2)+scale_colour_manual(values = my_palette) + theme(legend.position="none")
dev.off()



eigs = pca$sdev^2
var.percent = eigs / sum(eigs) * 100

pdf('var.percent.pdf')
barplot(var.percent, xlab="PC", ylab="Percent Variance", names.arg=1:length(var.percent), las=1, ylim=c(0,max(var.percent)), col="blue")
#abline(h=1/ncol(pca)*100, col="red")
dev.off()


pca_mat = pca$x

pdf('var.percent.0rep.pdf')
eigs_0rep = apply(pca_mat[c(1,12,29,30,31),], 2, function(x) sum(x^2)/(dim(pca_mat)[1]-1))
var.percent.rep0 = eigs_0rep / sum(eigs) * 100
barplot(var.percent.rep0, xlab="PC", ylab="Percent Variance", names.arg=1:length(var.percent), las=1, ylim=c(0,max(var.percent)), col="blue")
dev.off()


pdf('var.percent.rep1.pdf')
eigs_1rep = apply(pca_mat[c(2, 4, 6, 8, 10, 13, 15, 17, 19, 21, 23, 25, 27),], 2, function(x) sum(x^2)/(dim(pca_mat)[1]-1))
var.percent.rep1 = eigs_1rep / sum(eigs) * 100
barplot(var.percent.rep1, xlab="PC", ylab="Percent Variance", names.arg=1:length(var.percent), las=1, ylim=c(0,max(var.percent)), col="blue")
dev.off()

pdf('var.percent.rep2.pdf')
eigs_2rep = apply(pca_mat[c(3, 5, 7, 9, 11, 14, 16, 18, 20, 22, 24, 26, 28),], 2, function(x) sum(x^2)/(dim(pca_mat)[1]-1))
var.percent.rep2 = eigs_2rep / sum(eigs) * 100
barplot(var.percent.rep2, xlab="PC", ylab="Percent Variance", names.arg=1:length(var.percent), las=1, ylim=c(0,max(var.percent)), col="blue")
dev.off()


