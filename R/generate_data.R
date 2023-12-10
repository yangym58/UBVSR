#' Generate data
#'
#' @param p A number gives the probability whether a covariate has effect
#' @param h A number represents the heritability
#' @param n A scalar: number of observations
#' @param M A scalar: total number of covariates
#'
#' @return A list with the elements
#' \item{Y}{A n x 1 vector of outcomes}
#' \item{X}{A n x M matrix of the whole design matrix}
#' \item{beta}{A M x 1 vector for the true coefficients of each covariates}
#' \item{epsilon}{A n x 1 vector of random effects}
#' @export
#'
#' @examples generate_data(0.05, 0.5, 1000, 200)
generate_data <- function(p, h, n, M){
  X <- matrix(rnorm(n * M), nrow = n)
  X <- scale(X, scale = F)
  beta <- rbinom(M, 1, p) * rnorm(M)
  epsilon <- rnorm(n)
  lambda <- as.numeric(sqrt(h * var(epsilon) / ((1 - h) * var(X%*%beta))))
  Y <- lambda * (X %*% beta) + epsilon
  Y <- scale(Y, scale = F)
  return(list(Y = Y, X = X, beta = beta, epsilon = epsilon))
}
