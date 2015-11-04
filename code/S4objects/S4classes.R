
# S4 classes

#{{{ numbers class
setClass("leopard",contains="integer",
         representation(
             expected.survival.rate = "numeric",    # empirical survival estimates
             realised.survival.rate = "numeric",    # realised stochastic survival
             maternal.effect        = "matrix",     # multiplicative survival adjustments for maternal age
             litter.size            = "numeric",    # number of cubs born per female
             expected.birth.rate    = "numeric",    # expectations per female age category
             realised.birth.rate    = "numeric",    # realised stochastic birth rate per female age category
             realised.birth         = "matrix",     # realised births per maternal age class
             sex.ratio              = "numeric",    # fixed sex ratio
             names                  = "character"))

# initialisation method
setMethod("initialize","leopard",function(.Object) {
    
    .Object@.Data                  <- integer(14)
    .Object@expected.survival.rate <- numeric(14)
    .Object@realised.survival.rate <- numeric(14)
    .Object@maternal.effect        <- array(1, dim = c(2,5), dimnames = list(age.class = c("cub","juvenile"), maternal.age = c("f36","f48","f60","f72","f84")))
    .Object@litter.size            <- numeric(5)
    .Object@expected.birth.rate    <- numeric(5)
    .Object@realised.birth.rate    <- numeric(5)
    .Object@realised.birth         <- array(0, dim = c(2,5), dimnames = list(age.class = c("cub","juvenile"), maternal.age = c("f36","f48","f60","f72","f84")))
    .Object@sex.ratio              <- 0.5
    
    .Object@names <- c("nc","nj","saf","f36","f48","f60","f72","f84","sam","m36","m48","m60","m72","m84")
    
    .Object
})

# constructor function
leopard <- function(numbers = NULL, survival.rates = NULL, litter.sizes = NULL) {
    
    leopard.object <- new('leopard')
    
    if(!is.null(numbers)) 
        leopard.object@.Data[] <- as.integer(numbers)
    
    if(!is.null(survival.rates)) 
        leopard.object@expected.survival.rate[] <- survival.rates
    
    if(!is.null(litter.sizes)) 
        leopard.object@litter.size[] <- litter.sizes
      
    leopard.object
}

# assignment function
setMethod("[<-",
          signature(x = "leopard"),
          function (x, i, j, ..., value) 
          {
              x@.Data[i] <- as.integer(value)
              x
          }
)
#}}}

###

#{{{ quota class
setClass("quota",contains="integer",representation(
  hquota      = "numeric",
  cc          = "numeric",
  cc_star     = "numeric",
  sr          = "numeric",
  thr         = "numeric",
  kill        = "numeric",
  delta       = "numeric",
  days        = "numeric",
  amin        = "numeric",
  imp         = "numeric"))

setMethod("initialize","quota",function(.Object) {
  .Object@.Data    <- integer(1)
  .Object@hquota   <- vector('numeric',14)
  .Object@amin     <- 7
  .Object@cc       <- 0.001
  .Object@cc_star  <- .Object@cc
  .Object@days     <- vector('numeric',5)
  .Object@imp      <- vector('numeric',14)
  .Object
})
#}}}


#{{{ leopard population class
setClass("lpop",representation(n="array",k="array",quota="array",days="array"))

setMethod("initialize","lpop",function(.Object,nyr.proj,nreps) {
  
  .Object@n <- array(0,dim=c(14,nyr.proj,nreps))
  dimnames(.Object@n)[[1]] <- c('nc','nj','saf','f36','f48','f60','f72','f84','sam','m36','m48','m60','m72','m84')
  if(nyr.proj>1) dimnames(.Object@n)[[2]] <- seq(0,nyr.proj-1,1)
  if(nreps>1) dimnames(.Object@n)[[3]] <- 1:nreps
  
  .Object@k <- array(0,dim=c(14,nyr.proj,nreps))
  dimnames(.Object@k)[[1]] <- c('nc','nj','saf','f36','f48','f60','f72','f84','sam','m36','m48','m60','m72','m84')
  if(nyr.proj>1) dimnames(.Object@k)[[2]] <- seq(0,nyr.proj-1,1)
  if(nreps>1) dimnames(.Object@k)[[3]] <- 1:nreps
  
  .Object@quota <- array(0,dim=c(nyr.proj,nreps))
  if(nyr.proj>1) dimnames(.Object@quota)[[1]] <- seq(0,nyr.proj-1,1)
  if(nreps>1) dimnames(.Object@quota)[[2]] <- 1:nreps
  
  .Object@days <- array(0,dim=c(nyr.proj,nreps))
  if(nyr.proj>1) dimnames(.Object@days)[[1]] <- seq(0,nyr.proj-1,1)
  if(nreps>1) dimnames(.Object@days)[[2]] <- 1:nreps
  
  .Object
  
})
#}}}







