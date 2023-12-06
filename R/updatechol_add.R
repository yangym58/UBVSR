#' Update Cholesky decomposition when adding a covariate
#'
#' @param rtri A right triangular matrix
#' @param sigma2 A number: the variance of coefficient
#' @param X The design matrix
#' @param indvec A vector indicates covariate
#' @param newv A number indicates new covariate
#'
#' @return The updated right triangular matrix
#' @export
#'
#' @examples
#' rtri <- matrix(c(1,0,1,1), nrow = 2)
#' X <- matrix(c(1,0,1,1,1,1), nrow = 2)
#' updatechol_add(rtri, 1, X, c(1,2), 3)
updatechol_add<-function(rtri, sigma2, X, indvec, newv)
{ indvec <- c(indvec, newv)
nonzero <- length(indvec)
xtx <- crossprod(X[, indvec]) + diag(1 / sigma2, ncol = nonzero, nrow = nonzero)
ltri <- t(rtri)
y <- xtx[, nonzero]
b <- rep(0, nonzero)
b[1] <- y[1] / ltri[1, 1]
for (i in 2:(nonzero-1)) {
  tmp <- y[i]
  for(j in 1:(i-1)){
    tmp <- tmp - ltri[i, j] * b[j]
  }
  b[i] <- tmp / ltri[i, i]
}
b[nonzero] <- sqrt(y[nonzero]- sum(b ^ 2))
rtri <- cbind(rbind(rtri, 0), b)
names(rtri) <- NULL
return(rtri)
}





