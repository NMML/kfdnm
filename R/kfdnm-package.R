#' @title Fit integrated known-fate, dynamic N-mixture models
#' 
#' @description The \code{kfdnm} package provide some basic functions for fitting itegrated known-fate, dynamic N-mixture models. 
#' The functions can be used for both maximum-likelihood and Bayesian (via MCMC) inference. A Hidden Markov Model (HMM) formulation allows efficient sampling
#' of the abundance process in an MCMC or calculation of the likelihood for optimization.
#' 
#' 
#' \tabular{ll}{
#' Package: \tab crawl\cr 
#' Type: \tab Package\cr 
#' Version: \tab 1.901\cr 
#' Date: \tab October 23, 2014\cr 
#' License: \tab Unlimited \cr 
#' LazyLoad: \tab yes\cr }
#' 
#' @name kfdnm-package
#' @aliases kfdnm-package kfdnm
#' @docType package
#' @author Devin S. Johnson
#' 
#' Maintainer: Devin S. Johnson <devin.johnson@@noaa.gov>
#' @import Rcpp RcppArmadillo
#' @useDynLib kfdnm
NULL
