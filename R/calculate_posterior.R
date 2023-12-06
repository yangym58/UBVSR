#' Calculate posterior
#'
#' @param rtri A right triangular matrix
#' @param X Design matrix
#' @param Y A vector of outcome
#' @param s The vector of variance in each column of X
#' @param gam A vector indicates whether has effect
#' @param sigma2 The variance
#'
#' @return The posterior in log form
#' @export
#'
#' @examples
#' rtri <- matrix(c(1,0,0,1,1,0,1,1,1), ncol = 3)
#' X <- matrix(1:9, ncol = 3)
#' calculate_posterior(rtri, X, c(1,1,1), c(1,1,1), c(1,1,1), 1)
calculate_posterior <- function(rtri, X, Y, s, gam, sigma2)
{ M <- ncol(X)
n <- length(Y)
ind <- which(gam == 1)
v1 <- crossprod(X[, ind], Y)
v2 <- forwardsolve(t(rtri), v1)
betahat <- backsolve(rtri, v2)
aa <- 1 - (crossprod(Y, X[, ind]) %*% betahat) / (crossprod(Y) - n * mean(Y) ^ 2)
if(aa <= 0) {aa <- 0.01}
pygiven <- -sum(log(diag(rtri))) - (n / 2) * log(aa)
a <- sum(gam * s)
prior <- log(a / (1 + sigma2 * a) ^ 2)
return(pygiven + prior)

}
