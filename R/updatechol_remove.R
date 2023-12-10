#' Update Cholesky Decomposition when removing a covariate
#'
#' @param rtri A p x p right triangular matrix
#' @param indvec A p x 1 vector indicates the covariate
#' @param oldv A number indicates which covariate to be removed
#'
#' @return The new p x p triangular matrix
#' @export
#'
#' @examples
#' X <- matrix(rnorm(1000), ncol=20)
#' rtri <- chol(crossprod(X))
#' rtriprime <- updatechol_remove(rtri, 1:20, 5)
updatechol_remove <- function(rtri, indvec, oldv)
{  remove_ind <- which(indvec == oldv)
indvec <- indvec[-remove_ind]
rtri <- rtri[, -remove_ind, drop = F]
p <- nrow(rtri)
if(remove_ind == p){
  rtri <- rtri[-p, ]
} else{
  for(i in remove_ind:(p-1)){
    a <- rtri[i, i]
    b <- rtri[(i+1), i]
    if(b != 0){
      res <- given_rotation(a, b)
      G <- diag(1, nrow = p, ncol = p)
      G[i, i] <- G [(i + 1), (i + 1)] <- res$c
      G[i, (i + 1)] <- -res$s
      G[(i + 1), i] <- res$s
      rtri <- G %*% rtri
    }
  }
  rtri <- rtri[-p, ]
}
return(rtri)
}
