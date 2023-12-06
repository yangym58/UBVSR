#' Update Cholesky Decomposition when removing a covariate
#'
#' @param rtri A right triangular matrix
#' @param indvec A vector indicates the covariate
#' @param oldv A number indicates which covariate to be removed
#'
#' @return The new triangular matrix
#' @export
#'
#' @examples
#' rtri <- matrix(c(1,0,1,1), nrow = 2)
#' updatechol_remove(rtri, c(1,2), 1)
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
