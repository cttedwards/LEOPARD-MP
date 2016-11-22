loader <- function(name, path='C:/RESEARCH/LEOPARD-MP/results/') {
  do.call(load,list(file=paste(path,name,'.Rdata',sep=''),envir=globalenv()))
}
