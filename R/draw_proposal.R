#' Draw proposal in each iteration of MCMC
#'
#' @param m The number of trials
#' @param a A value satisfies a/m is the probability draw from the first case
#' @param gam A M x 1 vector indicates whether covariate has effect
#' @param sigma2 A value represents the variance
#' @param delta A value represents the variance in first case
#'
#' @return A list of elements
#' \item{sigma2prime}{A number: the proposal sigma2}
#' \item{gammaprime}{A M x 1 vector indicates the proposal gam}
#' \item{indic}{A value indicates the case of the proposal}
#' \item{indv}{A value indicates the index of the variate being changed}
#' @export
#'
#' @examples draw_proposal(5,1,c(1,0,1), 1, 1)
draw_proposal <- function(m, a, gam, sigma2, delta)
{  prop <- c((a / m), (1 - a / m) / 2, (1 - a / m) / 2)
indic <- which(rmultinom(1, 1, prop) == 1)
if(indic == 1){
  sigma2 <- exp(rnorm(1, mean = log(sigma2), sd = delta))
  indv <- 0
} else if(indic == 2){
  newind <- sample(which(gam == 0), 1)
  gam[newind] <- 1
  indv <- newind
} else {
  oldind <- sample(which(gam == 1), 1)
  gam[oldind] <- 0
  indv <- oldind
}
return(list(sigma2prime = sigma2, gammaprime = gam, indic = indic, indv = indv))
}
