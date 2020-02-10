args <- commandArgs(trailingOnly = TRUE)
intab <- args[1]

data=read.table(intab, sep="\t", header=FALSE, row.names=1)
cells=gsub("_.*", "", rownames(data))
cellOrder=c("LSK", "CMP", "GMP", "MEP", "G1E", "ER4", "CFUE", "ERY", "CFUMK", "iMK", "MON", "NEU", "NK", "TCD4", "TCD8", "B");
data=data[order(match(cells,cellOrder)),]
#insert row of NAs if want to split line 16 - 17
data=rbind(data[1:16,],NA,data[17:28,])
cells=cells[order(match(cells,cellOrder))]
cells=c(cells[1:16],NA, cells[17:28])
reps=rownames(data)
reps[17]=" "
cellColors=c("LSK"="#228B22", "CMP"="#EE7600", "GMP"="#008B8B", "MEP"="#CD6600", "CFUE"="#FF3030", "ERY"="#FF0000", "CFUMK"="#8B7355", "iMK"="#8B5A2B", "NEU"="#4F94CD", "MON"="#1874CD", "G1E"="#EE6363", "ER4"="#CD5555", "B"="#8B1C62", "TCD4"="#7A378B", "NK"="#8B008B", "TCD8"="#68228B")

library(ggplot2)
pdf("eryExp.pdf")
p=ggplot(data=data) + geom_line(aes(x=reps, y=V2, group=1)) + geom_line(aes(x=reps, y=V3, group=1)) + geom_line(aes(x=reps, y=V4, group=1))  + geom_point(stat="identity", aes(x=reps, y=V2, color=cells), size=2) + geom_point(stat="identity", aes(x=reps, y=V3, color=cells), size=2) + geom_point(stat="identity", aes(x=reps, y=V4, color=cells), size=2) + scale_x_discrete(limits=c(reps))+theme(axis.text.x=element_text(angle=90,hjust=1)) + scale_y_continuous(limits = c(0,400)) + scale_color_manual(values=ccol, name="Cell types") + labs(y='number of GO:1903706 "regulation of hemopoiesis" genes expressed', x='') 
p
dev.off()

q()

