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
heatmap.2(cor_matrix, dendrogram="row", col=my_palette, breaks=breaks, trace="none", cellnote=cr, notecol="black", revC=TRUE, notecex=.5, symm=TRUE, distfun=dist.pear, symkey=FALSE, margins = c(8,8))
dev.off()

pdf('hc_tree.pdf')
plot(hclust(as.dist(1-cor_matrix), method = "ward.D", members = NULL))
dev.off()



###### PCA analysis
library(ggplot2)

pca=prcomp(t(signal_matrix), center = FALSE, scale. = FALSE)
scores=data.frame(pca$x)

cols=ncol(scores)
cells=rownames(scores)#gsub("_.*", "", rownames(scores))
scores=data.frame(scores,cells)
colnames(scores)[cols+1]="cell"

my_palette = c("#8B1C62","#ff3030","#8b7355","#ee7600","#cd5555","#FF0000","#ee2c2c","#ee6363","#008b8b","#ff7f00","#228b22","#cd6600","#8b5a2b","#1874cd","#4f94cd","#8b008b","#7a378b","#68228b")

pdf('pca_pc1_pc2.pdf')
cell = colnames(signal_matrix)
ggplot(scores, aes(x=PC2, y=PC1))+geom_point(aes(color=cell), size=6)+geom_point(aes(color=cell), size=6)+geom_text(aes(label=rownames(scores)), size=3)+scale_colour_manual(values = my_palette)
dev.off()

pdf('pca_pc1_pc3.pdf')
cell = colnames(signal_matrix)
ggplot(scores, aes(x=PC3, y=PC1), col=my_palette)+geom_point(aes(color=cell), size=6)+geom_point(aes(color=cell), size=6)+geom_text(aes(label=rownames(scores)), size=3)+scale_colour_manual(values = my_palette)
dev.off()


pdf('pca_pc2_pc3.pdf')
cell = colnames(signal_matrix)
ggplot(scores, aes(x=PC2, y=PC3), col=my_palette)+geom_point(aes(color=cell), size=6)+geom_text(aes(label=rownames(scores)), size=3)+scale_colour_manual(values = my_palette)
dev.off()


eigs = pca$sdev^2
var.percent = eigs / sum(eigs) * 100

pdf('var.percent.pdf')
barplot(var.percent, xlab="PC", ylab="Percent Variance", names.arg=1:length(var.percent), las=1, ylim=c(0,max(var.percent)), col="blue")
abline(h=1/ncol(pca)*100, col="red")
dev.off()




