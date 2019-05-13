#' SRAT without adjustment for additional covariates
#' 
#' @param Y outcome variable
#' @param Z covariates for testing
#' @param n_ptb number of perturbations
#' 
#' @examples
#' # Data from wilcox.test
#' x <- c(0.80, 0.83, 1.89, 1.04, 1.45, 1.38, 1.91, 1.64, 0.73, 1.46)
#' y <- c(1.15, 0.88, 0.90, 0.74, 1.21)
#' wilcox.test(x, y)
#' # SRAT p-value
#' set.seed(123)
#' Y <- c(x, y)
#' Z <- rep(0:1, c(length(x), length(y)))
#' srat0(Y, Z)$p_value
#' 
#' @export
srat0 <- function(Y, Z, n_ptb = 1000) {
  Z <- as.matrix(Z)
  N <- NROW(Z)
  p <- NCOL(Z)
  stopifnot(length(Y) == N)
  S_ptb <- matrix(NA, n_ptb, p)
  Q_ptb <- rep(NA, n_ptb)
  V <- rep(1, N)
  alpha <- NULL
  alpha_ptb <- NULL
  e <- sum_I(Y, V / N) - sum_I(-Y, V / N)
  S <- c(t(Z) %*% (V * e) / N) 
  Q <- c(t(S) %*% S)
  for (i in 1:n_ptb) {
    V_ptb <- stats::rexp(N)
    e_ptb <- sum_I(Y, V_ptb / N) - sum_I(-Y, V_ptb / N)
    S_ptb[i, ] <- c(t(Z) %*% (V_ptb * e_ptb) / N)
    Q_ptb[i] <- c(t(S_ptb[i, ] - S) %*% (S_ptb[i, ] - S))
  }
  SIGMA <- stats::cov(S_ptb - matrix(S, n_ptb, length(S), byrow = TRUE))
  lambda <- pmax(eigen(SIGMA, symmetric = TRUE, only.values = TRUE)$values, 0)
  p_value <- kuonen(Q, lambda)$pval
  list(Y = Y, Z = Z, alpha = alpha, alpha_ptb = alpha_ptb, 
       S = S, S_ptb = S_ptb, Q = Q, Q_ptb = Q_ptb, 
       SIGMA = SIGMA, lambda = lambda, p_value = p_value)
}

