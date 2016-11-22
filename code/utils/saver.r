saver <- function(..., name, path='C:/RESEARCH/LEOPARD-MP/results/') {
  save(..., file=paste(path,name,'.Rdata',sep=''))
}
