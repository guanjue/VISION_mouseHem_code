library(UpSetR)

data = read.table('allpk_count_mat.txt', header = FALSE)

datset_list = read.table('bed_list_all.txt', header=FALSE)

#set.seed(2018)
#used_id = sample(dim(data)[1], 10000)
#data_s = as.data.frame(data[used_id,-1])

data_binary = data[,-1]
colnames(data_binary) = c('VISION','iChIP','SCR_FL','SCR_all')#as.matrix(datset_list[,2])

data_binary[data_binary>0] = 1
#data_binary_s = data_binary[,-c(8,9,10)]

#pdf('ccREs_intersect_num_all0.pdf', height=7, width=20)
#upset(data_binary_s[,], nsets = 10, nintersects = NA, mb.ratio = c(0.7, 0.3), text.scale=3, point.size = 5, group.by = 'degree', scale.intersections="identity", scale.sets="identity", order.by = c("freq", 'degree'),decreasing = c(TRUE, FALSE))
#dev.off()
pdf('ccREs_intersect_num_all0_GZX_AQW.pdf', height=5, width=8, onefile=FALSE)
upset(data_binary, nsets = 10, nintersects = NA, group.by = 'degree', scale.intersections="identity", scale.sets="identity", order.by = c("freq", 'degree'),decreasing = c(TRUE, FALSE), text.scale = c(1.3, 1.3, 1, 1, 1.5, 1.5), point.size = 3, line.size = 0.5)
dev.off()


