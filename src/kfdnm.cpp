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

arma::vec make_S_prob_vec(const int& from_N, const int& R, const int& to_N, const double& omega, 
const double& gamma, const int N_max){
  int m = std::min(from_N-R, to_N);
  arma::vec prob(m+1, fill::ones);
  for(int i=0; i<=m; i++){
    prob(i) = R::dbinom(i, from_N-R, omega, false) * R::dpois(to_N-i, gamma, false);
    }
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
Rcpp::List dnm_hmm(
  const Rcpp::IntegerVector& n, 
  const Rcpp::IntegerVector& R, 
  const Rcpp::IntegerVector& id, 
  const Rcpp::NumericVector& omega_dnm,  
  const Rcpp::NumericVector& gamma, 
  const Rcpp::NumericVector& p, 
  const int& N_max, 
  const bool& back_sample
  ){
    int I = n.size();
    arma::mat S_probs(I, N_max+1, fill::zeros);
    arma::mat phi(I, N_max+1, fill::zeros);
    arma::rowvec eta(N_max+1, fill::ones); 
    for(int j=1; j<=N_max; j++){eta(j) = 1/(1.0*j);}
    eta = arma::normalise(eta, 1);
    arma::cube Delta(N_max+1, N_max+1, I, fill::zeros);
    double ll_dnm = 0;
    Rcpp::LogicalVector miss_n = Rcpp::is_na(n);
    arma::mat P_mat(N_max+1, N_max+1);
    for(int i = 0; i<I; i++){
      if(miss_n[i]){
        P_mat.eye();
      } else{
        P_mat = make_P_mat(n[i], p[i], R[i], N_max);
      }
      if( i==0 || id[i]!=id[i-1] ){
        phi.row(i) = eta * P_mat;
        ll_dnm += log(arma::sum(phi.row(i)));
      } else{
        Delta.slice(i) = N_trans_mat(omega_dnm[i-1], gamma[i], R[i-1], N_max);
        phi.row(i) = phi.row(i-1) * Delta.slice(i) * P_mat;
        ll_dnm += log(arma::sum(phi.row(i)));
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
        if(id[j-1]!=id[j]){
          prob = phi.row(j-1).t();
        } else{
          prob = phi.row(j-1).t() % Delta.slice(j).col(N[j]);
        }
        N[j-1] = sample_du(prob) - 1;
        if(j>0 && id[j]==id[j-1]){
          S[j]=sample_du(make_S_prob_vec(N[j-1], R[j-1], N[j], omega_dnm[j-1], gamma[j], N_max))-1;
          G[j]=N[j]-S[j];
        }
      }
      return Rcpp::List::create(
        Rcpp::Named("logLik")=ll_dnm,
        Rcpp::Named("N")=N,
        Rcpp::Named("S")=S,
        Rcpp::Named("G")=G);
    } else{
      return Rcpp::List::create(
        Rcpp::Named("logLik")=ll_dnm
        );
    }    
  }
  
  // [[Rcpp::export]]
  Rcpp::List kf_hmm(
    const Rcpp::IntegerVector& Y,
    const Rcpp::NumericVector& omega,
    const Rcpp::IntegerVector& id)
    {
      int I = Y.size();
      double ll = 0;
      //Rcpp::NumericVector ll(I);
      arma::mat phi(I, 2, fill::zeros);
      phi(0,1)=1;
      Rcpp::LogicalVector miss_Y = Rcpp::is_na(Y);
      arma::mat P_mat(2, 2, fill::zeros);
      arma::mat Delta(2,2,fill::eye);
      for(int i=1; i<I; i++){
        if(id[i]!=id[i-1]){
          phi(i,1) = 1;
        } else{
          if(miss_Y[i]){
            P_mat.eye();
          } else {
            P_mat.zeros();
            P_mat(Y[i],Y[i]) = 1;
          }
          Delta(1,0) = 1-omega[i-1];
          Delta(1,1) = omega[i-1];
          phi.row(i) = phi.row(i-1) * Delta * P_mat;
          ll += log(arma::sum(phi.row(i)));
          //ll[i] = log(arma::sum(phi.row(i)));
          phi.row(i) = arma::normalise(phi.row(i),1);
        } // end of missing data ifelse
      } // end of loop
      return Rcpp::List::create(
        Rcpp::Named("logLik")=ll
        );
    }
    
    