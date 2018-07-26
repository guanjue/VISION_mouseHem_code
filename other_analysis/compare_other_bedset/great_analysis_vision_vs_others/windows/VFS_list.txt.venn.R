library(venneuler)
vc = venneuler(c('V'=41224, 'F'=6064, 'S'=261444, 'F&S'=10331, 'V&S'=63017, 'V&F'=10305, 'V&F&S'=90473))

vc$labels = c('', '', '')
vc$colors = c(1, 0.5, 0.2)

pdf('VFS_list.txt.pdf')
plot(vc)
dev.off()

