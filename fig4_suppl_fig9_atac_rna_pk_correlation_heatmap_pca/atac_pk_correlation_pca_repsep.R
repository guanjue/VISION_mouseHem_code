#setwd('/Users/universe/Documents/2018_BG/vision_mouse_analysis/ct_hc_pca')
setwd('/storage/home/gzx103/group/projects/vision/ct_hc_pca')

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
breaks = c(seq(.2,1,length=16))

dist.pear <- function(x) as.dist(1-x)
pdf('hc_heatmap.pdf')
#pheatmap(cor_matrix, color=my_colorbar, cluster_cols = TRUE,cluster_rows=TRUE,annotation_names_row=TRUE,annotation_names_col=TRUE,show_rownames=TRUE,show_colnames=TRUE,clustering_method = 'complete')
heatmap.2(cor_matrix, dendrogram="row", col=my_palette, breaks=breaks, trace="none", notecol="black", revC=TRUE, notecex=.5, symm=TRUE, distfun=dist.pear, symkey=FALSE, margins = c(8,8))
dev.off()

pca=prcomp(t(signal_matrix), center = FALSE, scale. = FALSE)
scores=data.frame(pca$x)
pk_loading = pca$rotation


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
abline(h=1/ncol(pca)*100, col="red")
dev.off()




