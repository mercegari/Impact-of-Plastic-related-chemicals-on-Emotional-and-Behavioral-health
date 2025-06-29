#--------------------------------------------------------#
#      				    MY LIST OF FUNCTIONS				           #
#--------------------------------------------------------#


##########################################################
# Funció CI
ci <- function(x, level=.95, log=FALSE) {
  if (log==TRUE) x <- log(x)
  mu <- mean(x)
  sigma <- sd(x)
  n <- length(x)
  error <- qnorm(1 - ((1-level)/2)) * sigma / sqrt(n)
  if (log==TRUE) {
    return(c(mean=exp(mu), ci.low=exp(mu-error), ci.high=exp(mu+error)))
  } else {
    return(c(mean=mu,ci.low=mu-error, ci.high=mu+error))
  }
}


# Escalar les variables a dues desviacions estàndards, com en Gelman:
escalamenta <- function(x) (x - mean(x, na.rm=TRUE))/(2*sd(x, na.rm=TRUE))


# Check percentage of missing data
# Extracted from 
# https://datascienceplus.com/imputing-missing-data-with-r-mice-package/

pMiss <- function(x){sum(is.na(x))/length(x)*100}

