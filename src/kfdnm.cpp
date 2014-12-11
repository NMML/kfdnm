// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

#include "RcppArmadillo.h"
// [[Rcpp::depends(RcppArmadillo)]]

using namespace arma;


//// [[Rcpp::export]]
//double N_trans_old(const int& j, const int& k, const double& omega, const double& gamma, const int& R){
//  double out=0;
//  if(j < R) return out;
//  for(int c=0; c<=std::min(j-R,k);  c++){
//    out += exp(R::dbinom(c, j-R, omega, true) + R::dpois(k-c, gamma, true));
//  }
//  return out;
//}

// [[Rcpp::export]]
arma::vec N_trans(const int& fromN, const double& omega, const double& gamma, 
const int& R, const int& N_max){
  Rcpp::IntegerVector N(N_max+1);
  for(int i=0; i<=N_max; i++){N(i)=i;}
  arma::vec a(Rcpp::dbinom(N, fromN-R, omega, 0));
  arma::vec b(Rcpp::dpois(N, gamma, 0));
  arma::vec c = conv(a,b);
  return(c.subvec(0,N_max));
}


// [[Rcpp::export]]
arma::mat N_trans_mat(const double& omega, const double& gamma, const int& R, const int& N_max) {
  arma::mat out(N_max+1, N_max+1, fill::zeros);
  for(int j=R; j<=N_max; j++){
    out.row(j) = N_trans(j,omega,gamma,R,N_max).t();
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

arma::vec make_S_prob_vec(const int& from_N, const int& R, const int& to_N, const double& omega, const double& gamma){
  Rcpp::IntegerVector S(from_N-R+1);
  for(int i=0; i<=from_N-R; i++){S[i]=i;}
  arma::vec prob(dbinom(S, from_N-R, omega, false) * dpois(to_N-S, gamma, false));
  return(prob);
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
  const Rcpp::IntegerVector& R, 
  const Rcpp::IntegerVector& new_group, 
  const Rcpp::NumericVector& omega_dnm, 
  const Rcpp::NumericVector& omega_kf, 
  const Rcpp::NumericVector& gamma, 
  const Rcpp::NumericVector& p, 
  const int& N_max, 
  const bool& back_sample){
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
        ll_kf += R::dbinom(Y[i], Y[i-1]+R[i-1], omega_kf[i], true);
      }
      phi.row(i) = arma::normalise(phi.row(i),1);
    }
    if(back_sample){
      Rcpp::IntegerVector N(I);
      Rcpp::IntegerVector S(I);
      std::fill( S.begin(), S.end(), Rcpp::NumericVector::get_na() ) ;
      Rcpp::IntegerVector G(I);
      std::fill(G.begin(), G.end(), Rcpp::NumericVector::get_na() ) ;
      arma::vec prob(N_max+1);
      N[I-1] = sample_du(phi.row(I-1).t());
      for(int j=I-1; j>0; j--){
        prob = phi.row(j-1).t() % Delta.slice(j).col(N[j]);
        N[j-1] = sample_du(prob) - 1;
        if(new_group[j]!=1){
          S[j]=sample_du(make_S_prob_vec(N[j-1], R[j-1], N[j], omega_dnm[j], gamma[j]))-1;
          G[j]=N[j]-S[j];
        }
      }
      return Rcpp::List::create(
        Rcpp::Named("n2ll_dnm")=-2*ll_dnm,
        Rcpp::Named("n2ll_kf")=-2*ll_kf,
        Rcpp::Named("N")=N,
        Rcpp::Named("S")=S,
        Rcpp::Named("G")=G);
    } else{
      return Rcpp::List::create(
        Rcpp::Named("n2ll_dnm")=-2*ll_dnm,
        Rcpp::Named("n2ll_kf")=-2*ll_kf);
    }   
    
  }
  
  
  