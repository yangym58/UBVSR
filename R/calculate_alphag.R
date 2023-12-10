#' Calculate alphag
#' @inheritParams calculate_posterior
#' @param indvec A p x 1 vector indicates covariate that have effects
#' @param post A value represents the former posterior
#' @param gammaprime A M x 1 represents the proposal gamma
#' @param sigma2prime A value represents proposal variance
#' @param indic A number indicates which case
#' @param indv A number indicates which variate is changed
#' @param g A number indicates which balance function
#'
#' @return A list with the elements
#' \item{rtriprime}{The p x p right triangular matrix for the propsal}
#' \item{alpha}{A value of the function alphag in log form}
#' @export
#'
#' @examples
#' rtri <- matrix(c(1,0,0,1,1,0,1,1,1), ncol = 3)
#' X <- matrix(1:9, ncol = 3)
#' calculate_alphag(rtri, X, c(1,1,1), c(1,1,1), c(1,2,3), 5, c(1,1,1), c(1,1,1), 1, 1.1, 1, 0, g = 1)
calculate_alphag <- function(rtri, X, Y, s, indvec, post, gam, gammaprime, sigma2, sigma2prime, indic, indv, g = 1){
  M <- ncol(X)
  if(g == 1){
    if(indic == 1){
      propratio <- sigma2prime / sigma2
      rtriprime <- chol(crossprod(X[, which(gammaprime == 1)]) + diag(1 / sigma2prime, nrow = sum(gammaprime), ncol = sum(gammaprime)))
      postprime <- calculate_posterior(rtriprime, X, Y, s, gammaprime, sigma2prime)
      alpha <- (log(propratio) + (-sum(gam) / 2) * log(sigma2prime / sigma2) + postprime - post) / 2
    } else if(indic == 2){
      propratio <- (M - sum(gam)) / (sum(gam) + 1)
      rtriprime <- updatechol_add(rtri, sigma2, X, indvec, indv)
      postprime <- calculate_posterior(rtriprime, X, Y, s, gammaprime, sigma2prime)
      alpha <- (log(sum(gam)) - log(M - sum(gam)) + log(propratio) - log(sqrt(sigma2)) + postprime - post) / 2
    } else{
      propratio <- sum(gam) / (M - sum(gam) + 1)
      rtriprime <- updatechol_remove(rtri, indvec, indv)
      postprime <- calculate_posterior(rtriprime, X, Y, s, gammaprime, sigma2prime)
      alpha <- (log(M + 1 -sum(gam)) - log(sum(gam) - 1) + log(propratio) + log(sqrt(sigma2)) + postprime - post) / 2
    }
  }
  return(list(rtriprime = rtriprime, alpha = alpha))
}
