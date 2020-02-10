library(pheatmap)

#setwd('/Users/universe/Documents/2018_BG/vision_mouse_analysis/ccREcounts_state_transition/')

ideas_state = read.table('atac_pk_no0.atac.ideas.txt', header=F)
state_mat = ideas_state[,c(-(12+4),-(16+4))]

state_mat_label = state_mat[,-c(1:4)]

ct_names = c('LSK', 'HPC7', 'CMP', 'MEP', 'G1E', 'ER4', 'CFUE', 'ERY', 'ERY_fl', 'CFUMK', 'iMK', 'GMP', 'MON', 'NEU', 'NK', 'B_SPL', 'T_CD4', 'T_CD8')

#for (a in c(1:18)){
for (a in c(1,3)){
#a = 1
#        for (b in c(1:18)){
        for (b in c(8,11)){


                if (a<b){
                cell_pair = paste(ct_names[a], ct_names[b], sep='_')
                print(cell_pair)
                ct1 = state_mat[,a+4]
                ct2 = state_mat[,b+4]

                state_pair = apply(cbind(ct1, ct2), 1, function(x) paste(x[1], x[2], sep='_'))

                ccRE_num = dim(state_mat_label)[1] * dim(state_mat_label)[2]
                ccRE_num_ct = dim(state_mat_label)[1] 

                state_pair_list = c()
                state_pair_exp = c()
                state_pair_obs = c()
                state_pair_enrich = c()
                k=0
                for (i in c(0:26)){
                        for (j in c(0:26)){
                                k = k+1
                                print(k)
                                state_pair_tmp = paste(i, j, sep='_')
                                state_pair_list[k] = state_pair_tmp
                                ### exp state pair count
                                mep_i = sum(state_mat_label==i)
                                mep_j = sum(state_mat_label==j)
                                exp = mep_i / ccRE_num * mep_j / ccRE_num * ccRE_num_ct
                                ### obs state pair count
                                obs = sum(state_pair==state_pair_tmp)
                                #print(exp)
                                state_pair_exp[k] = exp
                                state_pair_obs[k] = obs
                                state_pair_enrich[k] = (obs+100)/(exp+100)
                        }
                }

                ideas_state_pair = cbind(state_pair_list, state_pair_enrich, state_pair_obs, state_pair_exp)

                ideas_state_pair_mat = matrix(state_pair_list, nrow = 27, byrow = TRUE)

                ideas_state_pair_enrich_mat = matrix(state_pair_enrich, nrow = 27, byrow = TRUE)

                ideas_state_pair_obs_mat = matrix(state_pair_obs, nrow = 27, byrow = TRUE)

                colnames(ideas_state_pair_enrich_mat) = c(0:26)
                rownames(ideas_state_pair_enrich_mat) = c(0:26)

                colnames(ideas_state_pair_obs_mat) = c(0:26)
                rownames(ideas_state_pair_obs_mat) = c(0:26)


                #enrich_upperlim = 10
                #my_colorbar=colorRampPalette(c('white', 'red'))(n = 50)
                #col_breaks = c(seq(0, enrich_upperlim,length=33))
                #breaksList_enrich = seq(0, enrich_upperlim, by = 0.1)

                #png(paste('ideas_state_pair_enrich_mat.', cell_pair, '.png', sep=''))
                #pheatmap(ideas_state_pair_enrich_mat,cluster_rows = F,cluster_cols = F,color=my_colorbar,show_rownames=TRUE,show_colnames=TRUE,annotation_names_row=FALSE,annotation_names_col=TRUE, breaks = breaksList_enrich)
                #dev.off()

                enrich_upperlim = 5
                my_colorbar=colorRampPalette(c('white', 'red'))(n = 50)
                col_breaks = c(seq(0, enrich_upperlim,length=33))

                breaksList_enrich = seq(0, enrich_upperlim, by = 0.1)
                ideas_state_pair_enrich_mat0 = ideas_state_pair_enrich_mat
                #diag(ideas_state_pair_enrich_mat0) = 0
                ideas_state_pair_enrich_mat0[ideas_state_pair_enrich_mat0>enrich_upperlim]=enrich_upperlim
                #ideas_state_pair_enrich_mat0[ideas_state_pair_enrich_mat0>10] = 10
                pdf(paste('ideas_state_pair_enrich_mat.', cell_pair, '.modify.pdf', sep=''), width=4, height=4)
                pheatmap(ideas_state_pair_enrich_mat0,cluster_rows = F,cluster_cols = F,color=my_colorbar,show_rownames=TRUE,show_colnames=TRUE,annotation_names_row=FALSE,annotation_names_col=TRUE, breaks = breaksList_enrich)
                dev.off()

                #my_colorbar=colorRampPalette(c('white', 'red'))(n = 128)
                #col_breaks = c(seq(0, 25000,length=33))

                #breaksList_count = seq(0, 25000, by = 250)
                #diag(ideas_state_pair_obs_mat) = 1
                #ideas_state_pair_obs_mat[ideas_state_pair_obs_mat>25000]=25000
                #ideas_state_pair_obs_mat = log10(ideas_state_pair_obs_mat+1)
                #png(paste('ideas_state_pair_obs_mat.', cell_pair, '.png', sep=''))
                #pheatmap(ideas_state_pair_obs_mat,cluster_rows = F,cluster_cols = F,color=my_colorbar,show_rownames=TRUE,show_colnames=TRUE,annotation_names_row=FALSE,annotation_names_col=TRUE, breaks = breaksList_count)
                #dev.off()
                }
        }
}

