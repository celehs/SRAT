#' SRAT with adjustment for additional covariates
#' 
#' @param Y outcome vector (e.g., phenotype)
#' @param Z covariates for testing (e.g., genotype)
#' @param X covariates for adjustment (e.g., age, gender)
#' @param cluster family identifier
#' @param init initial parameters
#' @param w_sqrt square root of weight vector
#' @param bw bandwidth
#' @param bw_adj adjusted factor to bandwidth
#' @param n_ptb number of perturbations
#' 
#' @export
srat <- function(Y, Z, X, cluster = NULL, init = NULL, w_sqrt = NULL,
                  bw = NULL, bw_adj = NULL, n_ptb = 1000) {
  N <- NROW(Z)
  p <- NCOL(Z)
  D <- rep(1, N)
  stopifnot(length(Y) == N)
  if (is.null(cluster)) cluster <- 1:N
  id_tmp <- as.integer(factor(cluster))
  n <- length(unique(id_tmp))
  if (is.null(w_sqrt)) w_sqrt <- rep(1, p)
  S_ptb <- matrix(NA, n_ptb, p)
  Q_ptb <- rep(NA, n_ptb)
  V <- rep(1, N)
  X <- as.matrix(X)
  e0 <- sum_I(Y, V * D / n) - D * sum_I(-Y, V / n)
  if (is.null(init)) init <- stats::lm(e0 ~ X)$coef[-1] # no intercept
  init <- init/sum(abs(init))
  if (is.null(D)) {
    M <- obj_min(init, X, Y, V, V, TRUE, 1, 1e-6, 5000)
  } else {
    M <- obj_min(init, -X, -Y, D * V, V, TRUE, 1, 1e-6, 5000)
  }
  alpha <- M[1, -1]
  eta <- c(X %*% alpha)
  if (is.null(bw)) {
    bw_opt <- npregbw(rank(Y) ~ eta, bwmethod = "cv.aic")$bw
    # bw_opt <- npregbw(e0 ~ eta, bwmethod = "cv.aic")$bw
    if (is.null(bw_adj)) bw_adj <- 1
    bw <- bw_opt * bw_adj
  }
  e <- rank_res(Y, D, V, eta, bw)/n
  S <- c(t(Z) %*% (V * e) / n) * w_sqrt
  Q <- c(t(S) %*% S)
  alpha_ptb <- matrix(NA, n_ptb, length(alpha))
  for (i in 1:n_ptb) {
    V_ptb0 <- stats::rexp(n)
    V_ptb <- V_ptb0[id_tmp]
    M_ptb <- obj_min(alpha, X, Y, V_ptb, V_ptb, TRUE, 1, 1e-6, 5000)
    alpha_ptb[i, ] <- M_ptb[1, -1]
    eta_ptb <- c(X %*% alpha_ptb[i, ])
    e_ptb <- rank_res(Y, D, V_ptb, eta_ptb, bw)/n
    S_ptb[i, ] <- c(t(Z) %*% (V_ptb * e_ptb) / n) * w_sqrt
    Q_ptb[i] <- c(t(S_ptb[i, ] - S) %*% (S_ptb[i, ] - S))
  }
  SIGMA <- stats::cov(S_ptb - matrix(S, n_ptb, length(S), byrow = TRUE))
  lambda <- pmax(eigen(SIGMA, symmetric = TRUE, only.values = TRUE)$values, 0)
  pval <- kuonen(Q, lambda)$pval
  list(Y = Y, Z = Z, X = X, w_sqrt = w_sqrt, bw = bw, bw_adj = bw_adj,
       alpha = alpha, alpha_ptb = alpha_ptb, S = S, S_ptb = S_ptb,
       Q = Q, Q_ptb = Q_ptb, SIGMA = SIGMA, lambda = lambda, pval = pval)
}

#' SRAT with adjustment for additional covariates (null model)
#' 
#' @param Y outcome vector (e.g., phenotype)
#' @param X covariates for adjustment (e.g., age, gender)
#' @param cluster family identifier
#' @param init initial parameters
#' @param n_ptb number of perturbations
#' 
#' @export
skat.null <- function(Y, X, cluster = NULL, init = NULL, n_ptb = 1000) {
  N <- NROW(Z)
  p <- NCOL(Z)
  D <- rep(1, N)
  stopifnot(length(Y) == N)
  if (is.null(cluster)) cluster <- 1:N
  ID <- as.integer(factor(cluster))
  n <- length(unique(ID))
  V <- rep(1, N)
  X <- as.matrix(X)
  e0 <- sum_I(Y, V * D / n) - D * sum_I(-Y, V / n)
  if (is.null(init)) init <- stats::lm(e0 ~ X)$coef[-1] # no intercept
  init <- init/sum(abs(init))
  if (is.null(D)) {
    M <- obj_min(init, X, Y, V, V, TRUE, 1, 1e-6, 5000)
  } else {
    M <- obj_min(init, -X, -Y, D * V, V, TRUE, 1, 1e-6, 5000)
  }
  alpha <- M[1, -1]
  eta <- c(X %*% alpha)
  alpha_ptb <- matrix(NA, length(alpha), n_ptb)
  eta_ptb <- matrix(NA, N, n_ptb)
  V_ptb <- matrix(NA, N, n_ptb)
  for (i in 1:n_ptb) {
    V_ptb0 <- stats::rexp(n)
    V_ptb[, i] <- V_ptb0[ID]
    M_ptb <- obj_min(alpha, X, Y, V_ptb[, i], V_ptb[, i], TRUE, 1, 1e-6, 5000)
    alpha_ptb[, i] <- M_ptb[1, -1]
    eta_ptb[, i] <- c(X %*% alpha_ptb[i, ])
  }
  list(n = n, Y = Y, X = X, alpha = alpha, eta = eta, 
       alpha_ptb = alpha_ptb, eta_ptb = eta_ptb, V_ptb = V_ptb)  
}

#' Objective function minimization
#' 
#' @param coef initial parameter values
#' @param Z covariates matrix
#' @param y response vector
#' @param u a vector of weights
#' @param v a vector of weigths 
#' @param tree indicator of whether the tree algorithm should be used for fast ranking
#' @param nrep number of replications in the optimization algorithm
#' @param ftol function tolerance to determine convergence
#' @param nfun_max maximum number of function evaluations
#' 
#' @export
obj.min <- function(coef, Z, y, u = NULL, v = NULL, tree = TRUE, 
                    nrep = 1, ftol = 1e-6, nfun_max = 5000) {
  n <- length(y)
  if (is.null(u)) u <- rep(1, n)
  if (is.null(v)) v <- rep(1, n)
  M <- obj_min(coef, Z, y, u, v, tree, nrep, ftol, nfun_max)
  list(min = M[1, 1], param = M[1, -1])
}
