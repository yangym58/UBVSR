#' Calculate posterior
#'
#' @param rtri A p x p right triangular matrix
#' @param X A n x M Design matrix
#' @param Y A n x 1 vector of outcome
#' @param s A M x 1 vector for variance in each column of X
#' @param gam A M x 1 vector indicates whether covariate has effect
#' @param sigma2 A value represents the variance
#'
#' @return A value represents the posterior in log form
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
