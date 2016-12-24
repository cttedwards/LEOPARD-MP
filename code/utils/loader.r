loader <- function(name, path='../results/') {
  do.call(load,list(file=paste(path,name,'.Rdata',sep=''),envir=globalenv()))
}
