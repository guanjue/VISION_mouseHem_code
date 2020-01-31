args = commandArgs(trailingOnly=TRUE)

x = args[1]
y = args[2]
x_norm_output = args[3]

sx = scan(x)
sy = scan(y)

sx_n = (sx+1) / sum(sx) * sum(sy)

write.table(sx_n, x_norm_output, sep='\t', col.names=FALSE, row.names=FALSE, quote=FALSE)