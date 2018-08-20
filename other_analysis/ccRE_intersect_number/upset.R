library(UpSetR)

data = read.table('allpk_count_mat.txt', header = FALSE)

datset_list = read.table('bed_list_all.txt', header=FALSE)

#set.seed(2018)
#used_id = sample(dim(data)[1], 10000)
#data_s = as.data.frame(data[used_id,-1])

data_binary = data[,-1]
colnames(data_binary) = as.matrix(datset_list[,2])

data_binary[data_binary>0] = 1
#data_binary_s = data_binary[,-c(8,9,10)]
data_binary_s = data_binary[,]
print(head(data_binary_s))

pdf('ccREs_intersect_num_ep300.pdf', height=5, width=20)
#upset(data_binary_s[,c( 2,3,4, 6,7)], nsets = 5, nintersects = NA, group.by = 'degree', scale.intersections="identity", scale.sets="identity", order.by = c("freq", 'degree'),decreasing = c(TRUE, FALSE))
upset(data_binary_s[,c( 2,3,4,5,6)], nsets = 5, nintersects = NA, group.by = 'degree', scale.intersections="identity", scale.sets="identity", order.by = c("freq", 'degree'),decreasing = c(TRUE, FALSE))
dev.off()

pdf('ccREs_intersect_num_knownREs.pdf', height=5, width=20)
#upset(data_binary_s[,c(1,2,3,4, 6  )], nsets = 5, nintersects = NA, group.by = 'degree', scale.intersections="identity", scale.sets="identity", order.by = c("freq", 'degree'),decreasing = c(TRUE, FALSE))
upset(data_binary_s[,c(1,2,3,4,5 )], nsets = 5, nintersects = NA, group.by = 'degree', scale.intersections="identity", scale.sets="identity", order.by = c("freq", 'degree'),decreasing = c(TRUE, FALSE))
dev.off()


pdf('ccREs_intersect_num_all.pdf', height=5, width=20)
upset(data_binary_s[,], nsets = 10, nintersects = NA, group.by = 'degree', scale.intersections="identity", scale.sets="identity", order.by = c("freq", 'degree'),decreasing = c(TRUE, FALSE))
dev.off()

