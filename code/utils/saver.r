saver <- function(..., name, path='../results/') {
  save(..., file=paste(path,name,'.Rdata',sep=''))
}
