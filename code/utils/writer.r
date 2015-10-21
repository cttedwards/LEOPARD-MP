writer <- function(x, name, path='../data/', ...) {
  write.csv(x, file=paste(path,name,'.csv',sep=''),row.names=F, ...)
}
