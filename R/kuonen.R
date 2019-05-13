#' Kuonen's method to calculate p-value
#' 
#' @param q observed test statistic
#' @param lambda eigen-values of the empirical covariance matrix
#' @param tol tolorence level
#' @param maxit maximum number of iterations
#' @export
kuonen <- function(q, lambda, tol = 1e-8, maxit = 1000) {
  # The 'survey' package contains a similar function called 'saddle'
  K0 <- function(x) - 0.5 * sum(log(1 - 2 * x * lambda))
  K1 <- function(x) sum(lambda / (1 - 2 * x * lambda))
  K2 <- function(x) 2 * sum(lambda^2 / (1 - 2 * x * lambda)^2)
  lambda.max <- max(lambda)
  lambda <- lambda / lambda.max
  q <- q / lambda.max
  tmp <- hybrid(q, lambda, tol, maxit)
  niter <- tmp$niter
  x <- tmp$root
  if (niter == maxit) warning("Maximum number of iterations has been reached!")
  w <- sign(x) * sqrt(2 * (x * q - K0(x)))
  v <- x * sqrt(K2(x))
  pval <- stats::pnorm(w + log(v/w)/w, lower.tail = FALSE)
  list(pval = pval, niter = niter, fn.root = K1(x) - q)
}
