library(venneuler)
vc = venneuler(c('V'=75546, 'A'=9141, 'F'=13798, 'A&F'=1120, 'V&F'=72787, 'V&A'=28695, 'V&A&F'=27991))

vc$labels = c('', '', '')
vc$colors = c(1, 0.1, 0.5)

pdf('VAF_list.txt.pdf')
plot(vc)
dev.off()

