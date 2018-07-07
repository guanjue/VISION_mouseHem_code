euclidean_dist = read.table('IDEAS_state_meansignal_Euclidean_dist_matrix.txt', header=FALSE)

colnames(euclidean_dist) = c('B_SPL', 'CFU_E_ad', 'CFUMK', 'CLP', 'CMP', 'ER4', 'ERY_ad', 'ERY_fl', 'G1E', 'GMP', 'HPC7', 'LSK_BM', 'MEP', 'MK_imm_ad', 'MK_mat_fl', 'MONO_BM', 'NEU', 'NK_SPL', 'T_CD4_SPL', 'T_CD8_SPL')
rownames(euclidean_dist) = c('B_SPL', 'CFU_E_ad', 'CFUMK', 'CLP', 'CMP', 'ER4', 'ERY_ad', 'ERY_fl', 'G1E', 'GMP', 'HPC7', 'LSK_BM', 'MEP', 'MK_imm_ad', 'MK_mat_fl', 'MONO_BM', 'NEU', 'NK_SPL', 'T_CD4_SPL', 'T_CD8_SPL')

library(pheatmap)
my_colorbar=colorRampPalette(c('white', 'red'))(n = 128)
col_breaks = c(seq(0, 2000,length=33))


pdf('IDEAS_state_meansignal_Euclidean_dist_matrix.pdf')
pheatmap(euclidean_dist, color=my_colorbar, cluster_cols = TRUE,cluster_rows=TRUE,annotation_names_row=TRUE,annotation_names_col=TRUE,show_rownames=TRUE,show_colnames=TRUE,clustering_method = 'complete')
dev.off()


ed=round(euclidean_dist, digits=2)
library("gplots")
library("RColorBrewer")

dist.pear <- function(x) as.dist(x)
pdf('hc_heatmap.pdf')
#pheatmap(cor_matrix, color=my_colorbar, cluster_cols = TRUE,cluster_rows=TRUE,annotation_names_row=TRUE,annotation_names_col=TRUE,show_rownames=TRUE,show_colnames=TRUE,clustering_method = 'complete')
heatmap.2(euclidean_dist, dendrogram="row", col=my_colorbar, trace="none", cellnote="none", notecol="black", revC=TRUE, notecex=.5, symm=TRUE, distfun=dist.pear, symkey=FALSE, margins = c(8,8))
dev.off()

