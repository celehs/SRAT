#' Sequence Robust Association Test (SRAT)
#' 
#' SRAT (Dai, Yang, et al. 2017) <doi:10.1111/biom.12643> is a fully rank-based and 
#' flexible approach to test for association between a set of genetic variants and an outcome, 
#' while accounting for within-family correlation and adjusting for covariates. Comparing to 
#' existing methods, SRAT has the advantages of allowing for unknown correlation structures 
#' and weaker assumptions about the outcome distribution. 
#' 
#' @name SRAT-package
#' @aliases SRAT
#' @keywords package
#' @useDynLib SRAT
#' @import np Rcpp RcppArmadillo
NULL

