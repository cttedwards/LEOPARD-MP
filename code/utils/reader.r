reader <- function(name, path='../data/', ...) {
    read.csv(file=paste(path,name,'.csv',sep=''), ...)
}
