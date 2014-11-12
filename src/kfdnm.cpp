// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

#include "RcppArmadillo.h"
// [[Rcpp::depends(RcppArmadillo)]]

using namespace arma;


// [[Rcpp::export]]
double N_trans(const int& j, const int& k, const double& omega, const double& gamma, const int& R){
  double out=0;
  if(j < R) return out;
  for(int c=0; c<=std::min(j-R,k);  c++){
    out += exp(R::dbinom(c, j-R, omega, true) + R::dpois(k-c, gamma, true));
  }
  return out;
}


// [[Rcpp::export]]
arma::mat N_trans_mat(const double& omega, const double& gamma, const int& R, const int& N_max) {
  arma::mat out(N_max+1, N_max+1, fill::zeros);
  for(int j=R; j<=N_max; j++){
    for(int k=0; k<=N_max; k++){
      out(j,k) = N_trans(j,k,omega,gamma,R);
    }
  }
  return out;
}


arma::mat make_P_mat(const int& n, const double& p, const int& R, const int& N_max){
  arma::mat P_mat(N_max+1, N_max+1, fill::zeros);
  for(int j=R; j<=N_max; j++){
    P_mat(j,j) = R::dbinom(n,j-R,p,false);
  }
  return P_mat;
}


int sample_du(arma::vec ppp){
  arma::vec cdf = cumsum(arma::normalise(ppp,1));
  double U = Rcpp::as<double>(Rcpp::runif(1));
  int out = 1;
  if(U<= cdf[0]) return(out);
  else
  {
    for(int i=1; i<ppp.n_elem; i++){ 
      if(U <= cdf[i]){
        out = i+1;
        return(out);
      }
    }
    return(out);
  }
}

// [[Rcpp::export]]
Rcpp::List ikfdnm_hmm(
  const Rcpp::IntegerVector& n, 
  const Rcpp::IntegerVector& Y,
  const Rcpp::IntegerVector& M,
  const Rcpp::IntegerVector& R, 
  const Rcpp::IntegerVector& new_group, 
  const Rcpp::NumericVector& omega_dnm, 
  const Rcpp::NumericVector& omega_kf, 
  const Rcpp::NumericVector& gamma, 
  const Rcpp::NumericVector& p, 
  const int& N_max, 
  const bool& back_sample
  ){
    int I = n.size();
    arma::mat phi(I, N_max+1, fill::zeros);
    arma::rowvec eta(N_max+1, fill::ones); 
    for(int j=1; j<=N_max; j++){eta(j) = 1/(1.0*j);}
    eta = arma::normalise(eta, 1);
    arma::cube Delta(N_max+1, N_max+1, I, fill::zeros);
    double ll_dnm = 0;
    double ll_kf = 0;
    arma::mat P_mat(N_max+1, N_max+1);
    for(int i = 0; i<I; i++){
      P_mat = make_P_mat(n[i], p[i], R[i], N_max);
      if(new_group[i]==1){
        phi.row(i) = eta * P_mat;
        ll_dnm += log(arma::sum(phi.row(i)));
      } else{
        Delta.slice(i) = N_trans_mat(omega_dnm[i], gamma[i], R[i-1], N_max);
        phi.row(i) = phi.row(i-1) * Delta.slice(i) * P_mat;
        ll_dnm += log(arma::sum(phi.row(i)));
        ll_kf += R::dbinom(Y[i], M[i-1], omega_kf[i], true);
      }
      phi.row(i) = arma::normalise(phi.row(i),1);
    }
    arma::vec N(I, fill::zeros);
    if(back_sample){
      arma::vec prob(N_max+1);
      N[I-1] = sample_du(phi.row(I-1).t());
      for(int j=I-1; j>0; j--){
        prob = phi.row(j-1).t() % Delta.slice(j).col(N[j]);
        N[j-1] = sample_du(prob) - 1;
      }
    }   
    return Rcpp::List::create(
      Rcpp::Named("n2ll_dnm")=-2*ll_dnm,
      Rcpp::Named("n2ll_kf")=-2*ll_kf,
      Rcpp::Named("N")=N);
  }
  
  
  