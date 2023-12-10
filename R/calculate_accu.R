#' Calculate accuracy
#'
#' @param beta A M x 1 vector represents true vector of coefficients
#' @param gammamat A (niter+1) x M estimated gamma matrix
#'
#' @return A value represents estimated accuracy
#' @export
#'
#' @examples calculate_accu(c(1,1,0), matrix(c(1,1,1,1,0,0), nrow = 2))
calculate_accu <- function(beta, gammamat)
{ gammamat <- na.omit(gammamat)
return(mean(which(beta!=0) %in% which(gammamat[nrow(gammamat), ]!=0)))
}
