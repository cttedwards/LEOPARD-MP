
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
