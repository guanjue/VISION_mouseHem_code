
shuf -n 20000 mouse_allChr_wTAD_weRP_wAllComp_071819.txt > mouse_ccRE_r.txt


library(umap)

d = read.table('mouse_ccRE_r.txt', header=F, sep='\t')

eRPs = d[,c(c(10:30),c(31:70))]
eRPs = d[,c(41:48)]
eRPs[is.na(eRPs)] = 0

max(eRPs[,-c(1:20)])

ptm <- proc.time()
eRPs.umap = umap(eRPs)
proc.time() - ptm

correlations = as.numeric(as.character(d[,30]))
correlations[is.na(correlations)] = -0.1

rbPal <- colorRampPalette(c('blue', 'green', 'purple', 'orange', 'red'))
Col_points <- rbPal(100)[as.numeric(cut(correlations,breaks = 100))]

pdf('test.umap.pdf')
plot(eRPs.umap$layout[,1],eRPs.umap$layout[,2], pch=1, col=Col_points)
dev.off()



pdf('test.umap.pdf')
plot(eRPs.umap$layout[,1],eRPs.umap$layout[,2], pch=1, col='gray')
points(eRPs.umap$layout[d[,73]=='No',1],eRPs.umap$layout[d[,73]=='No',2], pch=1, col='blue')
points(eRPs.umap$layout[d[,73]=='Yes',1],eRPs.umap$layout[d[,73]=='Yes', 2], pch=1, col='red')
dev.off()


pdf('test.umap.pdf')
plot(eRPs.umap$layout[,1],eRPs.umap$layout[,2], pch=1, col='gray')
points(eRPs.umap$layout[d[,30]==1,1],eRPs.umap$layout[d[,30]==1,2], pch=1, col='blue')
points(eRPs.umap$layout[d[,30]==2,1],eRPs.umap$layout[d[,30]==2,2], pch=1, col='green')
points(eRPs.umap$layout[d[,30]==3,1],eRPs.umap$layout[d[,30]==3,2], pch=1, col='purple')
points(eRPs.umap$layout[d[,30]==4,1],eRPs.umap$layout[d[,30]==4,2], pch=1, col='orange')
points(eRPs.umap$layout[d[,30]==5,1],eRPs.umap$layout[d[,30]==5,2], pch=1, col='red')
dev.off()


chr	tss	geneID	ccRE_start	ccRE_end	correlation	selected	gene_group	pk	LSK_BM	HPC7	CMP	MEG1E	ER4	CFU_E_ad	ERY_ad	ERY_fl	CFUMK	MK_imm_ad	MK_mat_fl	GMP	MONO_BM	NEU	CLP	NK_SPL	B_SPL	T_CD4_SPL	T_CD8_SPL	level	B_SPL_tss	B_SPL_dist	CFU_E_ad_tss	CFU_E_ad_dist	CFUMK_tss	CFUMK_dist	CLP_tss	CLP_dist	CMP_tss	CMP_dist	ER4_tss	ER4_dist	ERY_ad_tss	ERY_ad_dist	ERY_fl_tss	ERY_fl_dist	G1E_tss	G1E_dist	GMP_tss	GMP_dist	HPC7_tss	HPC7_dist	LSK_BM_tss	LSK_BM_dist	MEP_tss	MEP_dist	MK_imm_ad_tss	MK_imm_ad_dist	MK_mat_fl_tss	MK_mat_fl_dist	MONO_BM_tss	MONO_BM_dist	NEU_tss	NEU_dist	NK_SPL_tss	NK_SPL_dist	T_CD4_SPL_tss	T_CD4_SPL_dist	T_CD8_SPL_tss	T_CD8_SPL_dist	G1E-ER4_ccRE_compartment	Gene_Compartment	SameCompartment