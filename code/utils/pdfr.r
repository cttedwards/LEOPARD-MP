pdfr <- function(x, ..., name, path='../report/') {
  pdf(..., file=paste(path,name,'.pdf',sep=''))
  suppressWarnings(print(x))
  dev.off()
}