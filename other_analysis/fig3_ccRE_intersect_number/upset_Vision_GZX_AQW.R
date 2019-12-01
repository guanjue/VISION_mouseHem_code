library(UpSetR)

### new input for Alex
data = read.table('allpk_count_mat_withcolnames.csv', header = TRUE)

### Old input and code used for reference
#datset_list = read.table('bed_list_all.txt', header=FALSE)

#set.seed(2018)
#used_id = sample(dim(data)[1], 10000)
#data_s = as.data.frame(data[used_id,-1])

#data_binary = data[,-1]
#colnames(data_binary) = as.matrix(datset_list[,2])

#data_binary[data_binary>0] = 1
#data_binary_s = data_binary[,-c(8,9,10)]
#data_binary_s = data_binary[,]
#print(head(data_binary_s))

#pdf('ccREs_intersect_num_all0.pdf', height=5, width=20)
#upset(data_binary_s[,], nsets = 10, nintersects = NA, group.by = 'degree', scale.intersections="identity", scale.sets="identity", order.by = c("freq", 'degree'),decreasing = c(TRUE, FALSE))
#dev.off()


### Important information #####
### 
### To change text sizes: (relative to 1)
### text.scale = c(intersection size title, intersection size tick labels, set size title, set size tick labels, set names, numbers above bars)
### 
### To change dot or line size: (relative to 1)
### point.size = 1, line.size = 1
###
### To change colors:
### shade.color = Color of row shading in matrix
### sets.bar.color = Color of set size bar plot
### main.bar.color = Color of the main bar plot
### matrix.color = Color of the intersection points
###
### For queries:
### queries = Unified querie of intersections, elements, and custom row functions. Entered as
###   a list that contains a list of queries. query is the type of query being conducted.
###   params are the parameters of the query (if any). color is the color of the points on
###   the plot that will represent the query. If no color is selected one will be provided
###   automatically. active takes TRUE or FALSE, and if TRUE, it will overlay the
###   bars present with the results from the query. If FALSE a tick mark will indicate
###   the intersection size. See examples section on how to do this.
### 
### query.legend = Position query legend on top or bottom of UpSet plot
###
### 

### To get rid of scientific notation: 
options(scipen=999)
###

data_binary = data[,-1]
data_binary[data_binary>0] = 1
colnames(data_binary) <- c("VISION", "iChIP", "SCR_FL", "SCR_all")
### Updated code to change readability
pdf('ccREs_intersect_num_all0_GZX_AQW.pdf', height=5, width=8, onefile=FALSE)
upset(data_binary, nsets = 10, nintersects = NA, group.by = 'degree', scale.intersections="identity", scale.sets="identity", order.by = c("freq", 'degree'),decreasing = c(TRUE, FALSE), text.scale = c(1.3, 1.3, 1, 1, 1.5, 1.5), point.size = 3, line.size = 0.5)
dev.off()
