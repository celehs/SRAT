// Author: Ming Yang (hustyangming@gmail.com)
// Reference: Abrevaya, J. Economics Letters 62 (1999) 279-285.

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
arma::mat bs_tree(const arma::vec & y) {
  vec uniy = unique(y); // return unique values sorted in ascending order
  int node, llim, rlim, i = 0, j = 0, m = uniy.n_elem;
  mat TREE(m, 3), NODE(m, 3, fill::zeros);
  TREE.fill(NA_REAL);
  TREE.col(0) = uniy;
  NODE(0, 0) = trunc((m - 1) / 2);
  NODE(0, 1) = 0;
  NODE(0, 2) = m - 1;
  while(i < m) {
    node = NODE(i, 0);
    llim = NODE(i, 1);
    rlim = NODE(i, 2);
    if(node > llim) {
      j++;
      NODE(j, 0) = trunc((llim + node - 1) / 2);
      NODE(j, 1) = llim;
      NODE(j, 2) = node - 1;
      TREE(node, 1) = NODE(j, 0);
    }
    if(node < rlim) {
      j++;
      NODE(j, 0) = trunc((node + 1 + rlim) / 2);
      NODE(j, 1) = node + 1;
      NODE(j, 2) = rlim;
      TREE(node, 2) = NODE(j, 0);
    }
    i++;
  }
  return TREE;
}

// [[Rcpp::export]]
double num_sum(const arma::vec & eta, const arma::vec & y,
               const arma::vec & u, const arma::vec & v, const arma::mat & TREE) {
  double val = 0.0;
  int i, node, n = eta.n_elem, m = TREE.n_rows;
  uvec idx = sort_index(eta);
  mat DATA(n, 10, fill::zeros);
  DATA.col(0) = eta(idx);
  DATA.col(1) = y(idx);
  DATA.col(2) = u(idx);
  DATA.col(3) = v(idx);
  for(i = 0; i < n; i++) {
    node = trunc((m - 1) / 2);
    while(DATA(i, 1) != TREE(node, 0)) {
      if(DATA(i, 1) < TREE(node, 0)) {
        DATA(node, 5) += DATA(i, 3);
        if(DATA(node, 7) == 0 || DATA(i, 0) > DATA(node, 6)) {
          DATA(node, 6) = DATA(i, 0);
          DATA(node, 7) = DATA(i, 3);
        } else {
          DATA(node, 7) += DATA(i, 3);
        }
        node = TREE(node, 1);
      } else {
        val += DATA(i, 2) * (DATA(node, 4)  + DATA(node, 5) -
          (DATA(i, 0) == DATA(node, 6)) * DATA(node, 7) -
          (DATA(i, 0) == DATA(node, 8)) * DATA(node, 9));
        node = TREE(node, 2);
      }
    }
    if(DATA(node, 9) == 0 || DATA(i, 0) > DATA(node, 8)) {
      DATA(node, 8) = DATA(i, 0);
      DATA(node, 9) = DATA(i, 3);
    } else {
      DATA(node, 9) += DATA(i, 3);
    }
    DATA(node, 4) += DATA(i, 3);
    val += DATA(i, 2) * (DATA(node, 5) -
      (DATA(i, 0) == DATA(node, 6)) * DATA(node, 7));
  }
  return val;
}

// [[Rcpp::export]]
double den_sum(const arma::vec & y, const arma::vec & u, const arma::vec & v) {
  double val = 0.0;
  int i, n = y.n_elem;
  uvec idx = sort_index(y);
  mat DATA(n, 5, fill::zeros);
  DATA.col(0) = y(idx);
  DATA.col(1) = u(idx);
  DATA.col(2) = v(idx);
  DATA.col(3) = DATA.col(2);
  for(i = 1; i < n; i++) {
    if(DATA(i, 0) == DATA(i - 1, 0)) {
      DATA(i, 3) += DATA(i - 1, 3);
      DATA(i, 4)  = DATA(i - 1, 4);
    } else {
      DATA(i, 4)  = DATA(i - 1, 4) + DATA(i - 1, 3);
    }
    val += DATA(i, 1) * DATA(i, 4);
  }
  return val;
}

// [[Rcpp::export]]
double obj_fun(const arma::rowvec & coef, const arma::mat & Z, const arma::vec & y,
               const arma::vec & u, const arma::vec & v, bool tree) {
  int i, j, n = y.n_elem;
  double num = 0.0, den = 0.0;
  mat TREE;
  vec eta = Z * coef.t();
  if(tree == true) {
    TREE = bs_tree(y);
    num += num_sum(eta, y, u, v, TREE);
    den += den_sum(y, u, v);
  } else {
    for(i = 0; i < n; i++) {
      if(u(i) != 0) {
        for(j = 0; j < n; j++) {
          if(y(i) > y(j)) {
            num += u(i) * v(j) * (eta(i) > eta(j));
            den += u(i) * v(j);
          }
        }
      }
    }
  }
  return -num / den;
}

// [[Rcpp::export]]
arma::rowvec rowvec_norm(const arma::rowvec & x, int p) {
  rowvec y;
  if(p == 0) {
    y = x;
  } else if(p == 1 || p == 2) {
    y = x / norm(x, p);
  } else {
    Rcpp::stop("Invalid 'p'!");
  }
  return y;
}

// [[Rcpp::export]]
double obj_try(arma::mat & M, const arma::uvec & idx, int ndim, int p, double fac, const arma::mat & Z,
               const arma::vec & y, const arma::vec & u, const arma::vec & v, bool tree) {
  double fac1 = (1 - fac) / ndim;
  double fac2 = fac1 - fac;
  rowvec rtry(ndim + 1);
  rtry(span(1, ndim)) = rowvec_norm(fac1 * sum(M.cols(1, ndim), 0) -
    fac2 * M(span(idx(ndim)), span(1, ndim)), p);
  rtry(0) = obj_fun(rtry(span(1, ndim)), Z, y, u, v, tree);
  if(rtry(0) < M(idx(ndim), 0)) {
    M.row(idx(ndim)) = rtry;
  }
  return rtry(0);
}

// [[Rcpp::export]]
arma::mat obj_min(const arma::rowvec & coef, const arma::mat & Z, const arma::vec & y,
                  const arma::vec & u, const arma::vec & v, bool tree, int nrep, double ftol, int nfun_max) {
  int i, j, nfun, ndim = coef.n_elem, p = 1;
  double ytmp, ytry, rtol;
  mat oldM, M(ndim + 1, ndim + 1, fill::zeros);
  rowvec rtry(ndim + 1);
  uvec idx(ndim + 1);
  M(span(0), span(1, ndim)) = rowvec_norm(coef, p);
  M(0, 0) = obj_fun(M(span(0), span(1, ndim)), Z, y, u, v, tree);
  for(i = 1; i < ndim + 1; i++) {
    M(i, i) = 1.0;
    M(span(i), span(1, ndim)) = rowvec_norm(coef + M(span(i), span(1, ndim)), p);
    M(i, 0) = obj_fun(M(span(i), span(1, ndim)), Z, y, u, v, tree);
  }
  for(j = 0; j < nrep; j++) {
    oldM = M;
    nfun = 0;
    while(1) {
      idx = sort_index(M.col(0));
      rtol = 2 * fabs(M(idx(ndim), 0) - M(idx(0), 0)) /
        (fabs(M(idx(ndim), 0)) + fabs(M(idx(0), 0)) + 1e-10);
      if(rtol < ftol) {
        M = M.rows(idx);
        break;
      }
      if(nfun > nfun_max) {
        Rcpp::stop("Maximum number of function evaluations exceeded!");
      }
      nfun += 2;
      ytry = obj_try(M, idx, ndim, p, -1.0, Z, y, u, v, tree);
      if(ytry < M(idx(0), 0)) {
        ytry = obj_try(M, idx, ndim, p, 2.0, Z, y, u, v, tree);
      } else if(ytry >= M(idx(ndim - 1), 0)) {
        ytmp = M(idx(ndim), 0);
        ytry = obj_try(M, idx, ndim, p, 0.5, Z, y, u, v, tree);
        if(ytry >= ytmp) {
          for(i = 0; i < ndim + 1; i++) {
            if((unsigned)i != idx(0)) {
              M(span(i), span(1, ndim)) = rowvec_norm(0.5 * M(span(i), span(1, ndim)) +
                0.5 * M(span(idx(0)), span(1, ndim)), p);
              M(i, 0) = obj_fun(M(span(i), span(1, ndim)), Z, y, u, v, tree);
            }
          }
          nfun += ndim;
        }
      } else {
        nfun -= 1;
      }
    }
    // Rcpp::Rcout << "j = " << j << ", nfun = " << nfun << std::endl;
    if(j < nrep - 1) {
      for(i = 0; i < ndim + 1; i++) {
        M(span(i), span(1, ndim)) = rowvec_norm(M(span(0), span(1, ndim)) +
          0.4 * (oldM(span(i), span(1, ndim)) -
          oldM(span(0), span(1, ndim))), p);
        M(i, 0) = obj_fun(M(span(i), span(1, ndim)), Z, y, u, v, tree);
      }
      ftol = 0.6 * ftol;
    }
  }
  M.col(0) = M.col(0);
  idx = sort_index(M.col(0));
  return M.rows(idx);
}

// [[Rcpp::export]]
arma::vec sum_I(const arma::vec & y, const arma::vec & w) {
  // objective: s[i] = sum_j I(y[i] > y[j]) * w[j]
  uvec idx = sort_index(y);
  int i, n = y.n_elem;
  double tmp = 0.0;
  vec s(n, fill::zeros);
  for (i = 1; i < n; i++) {
    tmp += w(idx(i - 1));
    if (y(idx(i)) > y(idx(i - 1))) {
      s(idx(i)) = s(idx(i - 1)) + tmp;
      tmp = 0.0;
    } else {
      s(idx(i)) = s(idx(i - 1));
    }
  }
  return s;
}

// [[Rcpp::export]]
arma::vec rank_res(const arma::vec & Y, const arma::vec & D,
                   const arma::vec & V, const arma::vec & eta, double h) {
  int i, j, N = Y.n_elem;
  vec res(N, fill::zeros);
  double c, u, K;
  c = sqrt(2.0 * datum::pi);
  for (i = 0; i < N; i++) {
    for (j = 0; j < N; j++) {
      if (Y(i) > Y(j)) {
        u = (eta(i) - eta(j)) / h;
        K = exp(- u * u / 2.0) / c / h;
        res(i) += D(j) * V(j) * K;
      } else if (Y(i) < Y(j)) {
        u = (eta(i) - eta(j)) / h;
        K = exp(- u * u / 2.0) / c / h;
        res(i) -= D(i) * V(j) * K;
      }
    }
  }
  return res;
}

// [[Rcpp::export]]
double f0(double x, const arma::vec & lambda, double q) {
  double ans = 0.0;
  int i, n = lambda.n_elem;
  for (i = 0; i < n; i++) {
    ans += lambda(i) / (1.0 - 2.0 * x * lambda(i));
  }
  ans -= q;
  return ans;
}

// [[Rcpp::export]]
double f1(double x, const arma::vec & lambda) {
  double tmp, ans = 0.0;
  int i, n = lambda.n_elem;
  for (i = 0; i < n; i++) {
    tmp = lambda(i) / (1.0 - 2.0 * x * lambda(i));
    ans += tmp * tmp;
  }
  ans *= 2.0;
  return ans;
}

// [[Rcpp::export]]
Rcpp::List hybrid(double q, const arma::vec & lambda, double tol, int maxit) {
  int i, niter;
  double tmp, err, x0, x1, a = -0.5, b = 0.5; // b = 0.5 (function value = Inf)
  while (f0(a, lambda, q) > 0.0) a -= b - a;
  // while (f0(b, lambda, q) < 0.0) b += b - a;
  x0 = (a + b) / 2.0;
  for (i = 0; i < maxit; i++) {
    x1 = x0 - f0(x0, lambda, q) / f1(x0, lambda);
    if ((x1 > a) & (x1 < b)) {
      err = fabs(x1 - x0) / (fabs(x1) + tol);
      if (err < tol) break;
      x0 = x1;
    } else {
      tmp = f0(x0, lambda, q);
      if (tmp < 0.0) a = x0;
      if (tmp > 0.0) b = x0;
      x0 = (a + b) / 2.0;
    }
  }
  niter = i + 1;
  return Rcpp::List::create(Rcpp::Named("niter") = niter,
                            Rcpp::Named("root") = x1);
}
