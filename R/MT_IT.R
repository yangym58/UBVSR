#' The MT_IT methods
#'
#' @param Y A vector of outcome
#' @param X The design matrix
#' @param m The number of trials
#' @param a a/m the probability in first case proposal
#' @param delta The variance in first case proposal
#' @param niter Number of iterations
#' @param g The balance function
#'
#' @return The MCMC output
#' @export
#'
#' @examples
#' data <- generate_data(0.05, 0.4, 1000, 200)
#'   Y <- data$Y
#'   X <- data$X
#'   m <- 20
#'   a <- 2
#'   delta <- 0.5
#'   niter <- 1500
#'   res <- MT_IT(Y, X, m, a, delta, niter)

MT_IT <- function(Y, X, m, a, delta, niter, g = 1){
  n <- length(Y)
  M <- ncol(X)
  s <- apply(X, 2, var)
  # initializing values?
  sigma2 <- 1
  gam <- rbinom(M, 1, prob = 0.1)
  indvec <- which(gam == 1)
  # save outputs
  gammamat <- matrix(nrow = (niter + 1), ncol = M)
  sigma2vec <- rep(0, (niter + 1))
  wmat <- rep(0, (niter + 1))
  gammamat[1, ] <- gam
  sigma2vec[1] <- sigma2
  # draw trials
  Ssigma2 <- rep(0, m)
  Sgam <- matrix(nrow = m, ncol = M)
  Sindic <- rep(0, m)
  Sindv <- rep(0, m)
  rtrilist <- list()
  alphag <- rep(0, m)
  for(j in 1:m){
    res <- draw_proposal(m, a, gam, sigma2, delta)
    Ssigma2[j] <- res$sigma2prime
    Sgam[j, ] <- res$gammaprime
    Sindic[j] <- res$indic
    Sindv[j] <- res$indv
  }
  # Obtain the first Cholesky decomposition
  rtri <- chol(crossprod(X[, which(gam == 1)]) + diag(1 / sigma2, nrow = sum(gam), ncol = sum(gam)))
  # iteration
  for(i in 1:niter){
    post <- calculate_posterior(rtri, X, Y, s, gam, sigma2)
    for(j in 1:m){
      res1 <- calculate_alphag(rtri, X, Y, s, indvec, post, gam, Sgam[j, ], sigma2, Ssigma2[j], Sindic[j], Sindv[j], g)
      alphag[j] <- res1$alpha
      rtrilist[[j]] <- res1$rtriprime
    }
    largestind <- which.max(alphag)
    z <- alphag[largestind] + log(sum(exp(alphag - alphag[largestind])))
    wmat[i] <- exp(-z)
    if(sum(is.na(exp(alphag - z))) != 0) {break}
    nextind <- which(rmultinom(1, 1, exp(alphag - z)) == 1)
    gamtemp <- gam
    sigma2temp <- sigma2
    gam <- gammamat[(i + 1), ] <- Sgam[nextind, ]
    if(sum(gam) == 1 || sum(gam) == ncol(X)){
      break
    }
    sigma2 <- sigma2vec[i + 1] <- Ssigma2[nextind]
    rtri <- rtrilist[[nextind]]
    if(Sindic[nextind] == 2){
      indvec <- c(indvec, Sindv[nextind])
      Sindic[1] <- 3
      Sindv[1] <- Sindv[nextind]
    } else if(Sindic[nextind] == 3){
      indvec <- indvec[-which(indvec == Sindv[nextind])]
      Sindic[1] <- 2
      Sindv[1] <- Sindv[nextind]
    } else{
      Sindic[1] <- 1
      Sindv[1] <- 0
    }
    Ssigma2[1] <- sigma2temp
    Sgam[1, ] <- gamtemp
    for(j in 2:m){
      res <- draw_proposal(m, a, gam, sigma2, delta)
      Ssigma2[j] <- res$sigma2prime
      Sgam[j, ] <- res$gammaprime
      Sindic[j] <- res$indic
      Sindv[j] <- res$indv
    }
  }
  return(list(gammamat = gammamat, sigma2vec = sigma2vec, wmat = wmat))
}
