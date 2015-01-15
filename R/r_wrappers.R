#' @title Create IKFDNM model likehood function
#' 
#' @description This function assimilates the data and desirted model to create a likelihood function that is a function with only one argument, \code{par},
#' This likelihood function can teither be optimized via \code{\link{optim}} or placed within an MCMC routine.
#' 
#' @param survival A model statement for the survival parameters of the N-mixture portion of the model
#' @param recruit A model statement for the recruitment portion of the N-mixture model
#' @param detection A model statement for the the detection portion of the N-mixture model
#' @param kf_survival_effects A model statement describing the difference between the known-fate survival and the survival of the N-mixture described in 
#' \code{dnm_survival}. 
#' @param fixed_list A list object with named entries: \code{survival}, \code{recruit}, and \code{det}. Each of the entries is a vector of length 
#' \code{nrow(data)}.
#' @param data A data.frame object containing the data for analysis.
#' @param N_max An integer giving the maximum abundance of the N-mixture population.
#' 
#' @details The data provided in the \code{data} agument contains rows for each survey within each group. The rows must be forst ordered by group, 
#' then time within group. The data must contain the following \code{colnames}: 
#' \itemize{
#' \item \code{group} A factor variable indicating the survey group.
#' \item \code{n} The observed abundance count at each survey time.
#' \item \code{R} The number of new KF individuals sampled from the population at each time.
#' }
#' 
#' The \code{kf_survival_effects} model state provides a method for specifying the \emph{difference} in survival between the N-mixtue individuals and the 
#' known-fate individuals. Thus, the known-fate individual survival model is equal to \code{~survival_dnm + kf_survival_effects}.
#'     
#' @return A function of the vector \code{par}
#' @author Devin S. Johnson
#' @export
make_ikfdnm_likelihood = function(survival=~1, recruit=~1, detection=~1, kf_survival_effects=NULL, 
                                  fixed_list=NULL, kf_data, dnm_data, N_max, ln_prior=NULL){
  X_dnm = model.matrix(survival, dnm_data)
  X_kf = model.matrix(survival, kf_data)
  if(!is.null(kf_survival_effects)){
    X_kf = cbind(X_kf,model.matrix(kf_survival_effects, kf_data))
  }
  if(!is.null(fixed_list$kf_survival)) {X_kf[!is.na(fixed_list$kf_survival),]=0}
  if(!is.null(fixed_list$dnm_survival)) {X_dnm[!is.na(fixed_list$dnm_survival),]=0}
  X_kf = as.matrix(X_kf[,colSums(X_kf)!=0])
  X_dnm = as.matrix(X_dnm[,colSums(X_dnm)!=0])
  W = model.matrix(recruit, dnm_data)
  if(!is.null(fixed_list$recruit)){W[!is.na(fixed_list$recruit),]=0}
  W = as.matrix(W[,colSums(W)!=0])
  Z = model.matrix(detection, dnm_data)
  if(!is.null(fixed_list$detection)){Z[!is.na(fixed_list$detection),]=0}
  Z = as.matrix(Z[,colSums(Z)!=0])
  np_dnm = ncol(X_dnm)
  np_kf = ncol(X_kf)
  np_rec = ncol(W)
  np_det = ncol(Z)
  n2ll = function(par){
    beta_kf=par[1:np_kf]
    beta_dnm = par[1:np_dnm]
    rho=par[(np_kf+1):(np_kf+np_rec)]
    alpha=par[(np_kf+np_rec+1):(np_kf+np_rec+np_det)]
    omega_dnm = plogis(X_dnm%*%beta_dnm)
    if(!is.null(fixed_list$dnm_survival)){
      omega_dnm=ifelse(is.na(fixed_list$dnm_survival), omega_dnm, fixed_list$dnm_survival)
    }
    omega_kf = plogis(X_kf%*%beta_kf)
    if(!is.null(fixed_list$kf_survival)){
      omega_kf=ifelse(is.na(fixed_list$kf_survival), omega_kf, fixed_list$kf_survival)
    }
    gamma = exp(W%*%rho)
    if(!is.null(fixed_list$recruit)){
      gamma=ifelse(is.na(fixed_list$recruit), gamma, fixed_list$recruit)
    }
    p=plogis(Z%*%alpha)
    if(!is.null(fixed_list$detection)){
      p=ifelse(is.na(fixed_list$detection), p, fixed_list$detection)
    }
    R=ifelse(is.na(dnm_data$R), 0, dnm_data$R)
    ll_dnm = kfdnm:::dnm_hmm(
      n=dnm_data$n, 
      R=R,
      id=as.numeric(dnm_data$group), 
      omega_dnm=omega_dnm, 
      gamma=gamma,
      p=p, 
      N_max=N_max, 
      back_sample=FALSE)$logLik
    ll_kf = kfdnm:::kf_hmm(
      Y = kf_data$CH,
      omega=omega_kf,
      id=as.numeric(kf_data$id))$logLik
    if(!is.null(ln_prior)){
      ln_p = ln_prior(par)
    } else {
      ln_p = 0
    }
    return(-2*(ll_dnm+ll_kf+ln_p))
  }
  attr(n2ll,"npar") = np_kf + np_rec+np_det
  return(n2ll)
}

#' @title Create IKFDNM model function to sample recruits and survivors
#' 
#' @description This function assimilates the data and desirted model to create a function that will sample the number
#' of recruits and survivors given a parameter vector \code{par}.
#' 
#' @param survival A model statement for the survival parameters of the N-mixture portion of the model
#' @param recruit A model statement for the recruitment portion of the N-mixture model
#' @param detection A model statement for the the detection portion of the N-mixture model
#' @param kf_survival_effects A model statement describing the difference between the known-fate survival and the survival of the N-mixture described in 
#' \code{dnm_survival}. 
#' @param fixed_list A list object with named entries: \code{survival}, \code{recruit}, and \code{det}. Each of the entries is a vector of length 
#' \code{nrow(data)}.
#' @param data A data.frame object containing the data for analysis.
#' @param N_max An integer giving the maximum abundance of the N-mixture population.
#' 
#' @details The data provided in the \code{data} agument contains rows for each survey within each group. The rows must be forst ordered by group, 
#' then time within group. The data must contain the following \code{colnames}: 
#' \itemize{
#' \item \code{group} A factor variable indicating the survey group.
#' \item \code{n} The observed abundance count at each survey time.
#' \item \code{R} The number of new KF individuals sampled from the population at each time.
#' }
#' 
#' The \code{kf_survival_effects} model state provides a method for specifying the \emph{difference} in survival between the N-mixtue individuals and the 
#' known-fate individuals. Thus, the known-fate individual survival model is equal to \code{~survival_dnm + kf_survival_effects}.
#'     
#' @return A function of the vector \code{par}
#' @author Devin S. Johnson
#' @export
make_ikfdnm_sampler = function(survival=~1, recruit=~1, detection=~1, kf_survival_effects=NULL, 
                                  fixed_list=NULL, kf_data, dnm_data, N_max, ln_prior=NULL){
  X_dnm = model.matrix(survival, dnm_data)
  X_kf = model.matrix(survival, kf_data)
  if(!is.null(kf_survival_effects)){
    X_kf = cbind(X_kf,model.matrix(kf_survival_effects, kf_data))
  }
  if(!is.null(fixed_list$kf_survival)) {X_kf[!is.na(fixed_list$kf_survival),]=0}
  if(!is.null(fixed_list$dnm_survival)) {X_dnm[!is.na(fixed_list$dnm_survival),]=0}
  X_kf = as.matrix(X_kf[,colSums(X_kf)!=0])
  X_dnm = as.matrix(X_dnm[,colSums(X_dnm)!=0])
  W = model.matrix(recruit, dnm_data)
  if(!is.null(fixed_list$recruit)){W[!is.na(fixed_list$recruit),]=0}
  W = as.matrix(W[,colSums(W)!=0])
  Z = model.matrix(detection, dnm_data)
  if(!is.null(fixed_list$detection)){Z[!is.na(fixed_list$detection),]=0}
  Z = as.matrix(Z[,colSums(Z)!=0])
  np_dnm = ncol(X_dnm)
  np_kf = ncol(X_kf)
  np_rec = ncol(W)
  np_det = ncol(Z)
 sampler = function(par){
    beta_kf=par[1:np_kf]
    beta_dnm = par[1:np_dnm]
    rho=par[(np_kf+1):(np_kf+np_rec)]
    alpha=par[(np_kf+np_rec+1):(np_kf+np_rec+np_det)]
    omega_dnm = plogis(X_dnm%*%beta_dnm)
    if(!is.null(fixed_list$dnm_survival)){
      omega_dnm=ifelse(is.na(fixed_list$dnm_survival), omega_dnm, fixed_list$dnm_survival)
    }
    omega_kf = plogis(X_kf%*%beta_kf)
    if(!is.null(fixed_list$kf_survival)){
      omega_kf=ifelse(is.na(fixed_list$kf_survival), omega_kf, fixed_list$kf_survival)
    }
    gamma = exp(W%*%rho)
    if(!is.null(fixed_list$recruit)){
      gamma=ifelse(is.na(fixed_list$recruit), gamma, fixed_list$recruit)
    }
    p=plogis(Z%*%alpha)
    if(!is.null(fixed_list$detection)){
      p=ifelse(is.na(fixed_list$detection), p, fixed_list$detection)
    }
    R=ifelse(is.na(dnm_data$R), 0, dnm_data$R)
    smp_dnm = kfdnm:::dnm_hmm(
      n=dnm_data$n, 
      R=R,
      id=as.numeric(dnm_data$group), 
      omega_dnm=omega_dnm, 
      gamma=gamma,
      p=p, 
      N_max=N_max, 
      back_sample=TRUE)
    ll_kf = kfdnm:::kf_hmm(
      Y = kf_data$CH,
      omega=omega_kf,
      id=as.numeric(kf_data$id))$logLik
    if(!is.null(ln_prior)){
      ln_p = ln_prior(par)
    } else {
      ln_p = 0
    }
    return(list(smp_dnm=smp_dnm, ll_kf=ll_kf))
  }
  attr(sampler,"npar") = np_kf + np_rec+np_det
  return(sampler)
}
