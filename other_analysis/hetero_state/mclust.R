library(mclust)
library(pheatmap)

d = read.table('G1E_state_27_signal.txt', header=F)

colnames(d) = c('ATAC', 'CTCF', 'H3K27ac', 'H3K27me3', 'H3K36me3', 'H3K4me1', 'H3K4me3', 'H3K9me3')

set.seed(2018)
BIC = mclustBIC(d)
print(BIC)
pdf('BIC.pdf')
plot(BIC)
dev.off()

set.seed(2018)
fit = Mclust(d, x = BIC)

nr1= dim(table(fit$classification))
fit_mean = c()
for (i in c(1:nr1)){
	sig_tmp = d[fit$classification==i,]
	sig_tmp_mean = mean(as.matrix(sig_tmp))
	fit_mean[i] = sig_tmp_mean
}
plot_cluster_rank = order(fit_mean)
plot_cluster_rank

dr_mclust_plot = c()
dr_mclust_plot_label = c()
for (i in plot_cluster_rank){
	print(sum(fit$classification==i))
	dr_mclust_plot = rbind(dr_mclust_plot, d[fit$classification==i,])
	dr_mclust_plot_label = rbind(dr_mclust_plot_label, cbind(rep(i, sum(fit$classification==i)), rep(i, sum(fit$classification==i))) )
}

my_colorbar=colorRampPalette(c('white', 'red'))(n = 128)

pdf('G1E.state.27.pdf')
pheatmap((dr_mclust_plot), color=my_colorbar, cluster_rows = FALSE, cluster_cols = FALSE,annotation_names_row=FALSE,annotation_names_col=FALSE,show_rownames=FALSE,show_colnames=FALSE)
dev.off()

pdf('G1E.state.27.label.pdf')
pheatmap(dr_mclust_plot_label, cluster_rows = FALSE, cluster_cols = FALSE,annotation_names_row=FALSE,annotation_names_col=FALSE,show_rownames=FALSE,show_colnames=FALSE)
dev.off()


[1] 370
[1] 219
[1] 298
[1] 375
[1] 368
[1] 297
[1] 281



