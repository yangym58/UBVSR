#' Draw proposal in each iteration of MCMC
#'
#' @param m The number of trials
#' @param a a/m is the probability draw from the first case
#' @param gam A vector indicates whether has effect
#' @param sigma2 The variance
#' @param delta The variance in first case
#'
#' @return A propsal
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
