#' Generate data
#'
#' @param p The probability whether a covariate has effect
#' @param h The heritability
#' @param n Number of observations
#' @param M Total number of covariate
#'
#' @return Y, X, a vector of coefficients, random effect
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
