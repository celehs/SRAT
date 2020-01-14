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
#' @param y outcome vector (e.g., phenotype)
#' @param X covariates for adjustment (e.g., age, gender)
#' @param B number of perturbations
#' 
#' @export
srat.null <- function(y, X, B = 1000) {
  X <- as.matrix(X)  
  n <- length(y)
  d <- rep(1, n)  
  v <- rep(1, n)
  e0 <- sum_I(y, v / n) - sum_I(-y, v / n)
  coef <- stats::lm(e0 ~ X)$coef[-1]
  init <- coef / sum(abs(coef))
  alpha <- obj_min(init, X, y, v, v, 
                   TRUE, 1, 1e-6, 5000)[1, -1]
  eta <- c(X %*% alpha)
  bw <- npregbw(rank(y) ~ eta, bwmethod = "cv.aic")$bw
  e <- c(rank_res(y, d, v, eta, bw) / n)
  alpha.ptb <- matrix(NA, length(alpha), B)
  v.ptb <- matrix(NA, n, B)
  e.ptb <- matrix(NA, n, B)
  for (b in 1:B) {
    v.ptb[, b] <- stats::rexp(n)
    alpha.ptb[, b] <- obj_min(alpha, X, y, v.ptb[, b], v.ptb[, b], 
                              TRUE, 1, 1e-6, 5000)[1, -1]
    e.ptb[, b] <- c(rank_res(y, d, v.ptb[, b], c(X %*% alpha.ptb[, b]), bw) / n)    
  }
  list(e = e, e.ptb = e.ptb, v.ptb = v.ptb)
}

#' SRAT with adjustment for additional covariates (testing)
#' 
#' @param obj object from null model
#' @param Z covariates for testing (e.g., genotype)
#' @param w.sqrt square root of weight vector
#' 
#' @export
srat.test <- function(obj, Z, w.sqrt = NULL) { 
  if (is.null(w.sqrt)) w.sqrt <- rep(1.0, NCOL(Z))  
  n <- NROW(Z)
  B <- NCOL(obj$v.ptb)
  score <- c(w.sqrt * t(Z) %*% obj$e) / n 
  score.ptb <- w.sqrt * t(Z) %*% (obj$v.ptb * obj$e.ptb) / n  
  SIGMA <- stats::cov(t(score.ptb - matrix(score, length(score), B)))
  lambda <- pmax(eigen(SIGMA, symmetric = TRUE, only.values = TRUE)$values, 0)
  kuonen(q = c(t(score) %*% score), lambda)$pval  
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
