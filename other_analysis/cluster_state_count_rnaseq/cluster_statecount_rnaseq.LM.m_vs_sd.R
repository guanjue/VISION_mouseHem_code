m = apply(gene_exp, 1, mean)
s = apply(gene_exp, 1, sd)


library(mclust)

m_vs_sd = (cbind(m, s))
BIC <- mclustBIC(m_vs_sd)
summary(BIC)
mod1 <- Mclust(scale(m_vs_sd), G = 4)
summary(mod1, parameters = TRUE)
c = mod1$classification

png('test.m_vs_sd.png', width=1000, height=1000)
par(mfrow=c(2,2))
plot(m, s)
abline(v=-4, col='red')
abline(h=2, col='red')
plot(m, s)
points(m[c==1], s[c==1], col='red')
points(m[c==2], s[c==2],col='blue')
points(m[c==3], s[c==3],col='green')
points(m[c==4], s[c==4],col='gray')
abline(v=-4, col='red')
abline(h=2, col='red')
hist(m, breaks=50)
abline(v=-4, col='red')
box()
hist(s, breaks=50)
abline(v=2, col='red')
box()
dev.off()


dr0_all
gene_exp_genegroup_vec
used_id_RF = sample(dim(dr0_all)[1], 20000)
modelRF = randomForest(dr0_all[used_id_RF,-55], as.factor(gene_exp_genegroup_vec[used_id_RF]), ntree = 500, mtry = 6, importance = TRUE)

fit_RF_pre = predict(modelRF, dr0_all[-used_id_RF,-55])

RF_ARI = adjustedRandIndex(fit_RF_pre, as.factor(gene_exp_genegroup_vec[-used_id_RF]))
RF_ARI


sum((fit_RF_pre==as.factor(gene_exp_genegroup_vec[-used_id_RF])))

table((fit_RF_pre, as.factor(gene_exp_genegroup_vec[-used_id_RF])))

fit_RF_R2pre = R2(fit_RF_pre, dr0_all_testing[used_id_chr3_4,55])


fit_RF_pre = predict(modelRF, dr0_all_testing[used_id_chr3_4,-55])
fit_RF_R2pre = R2(fit_RF_pre, dr0_all_testing[used_id_chr3_4,55])
fit_RF_R2pre



