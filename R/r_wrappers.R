#' @title Create IKFDNM model likehood function
#' 
#' @description This function assimilates the data and desirted model to create a likelihood function that is a function with only one argument, \code{par},
#' This likelihood function can teither be optimized via \code{\link{optim}} or placed within an MCMC routine.
#' 
#' @param dnm_survival A model statement for the survival parameters of the N-mixture portion of the model
#' @param dnm_recruit A model statement for the recruitment portion of the N-mixture model
#' @param dnm_det A model statement for the the detection portion of the N-mixture model
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
#' \item \code{num_released} The number of known-fate individuals released into the population at the time of the survey
#' \item \code{num_returns} The number of known-fate individuals observed at each survey time. These are a binomial sample of the \code{num_release} of the 
#' previous row.
#' \item \code{n} The observed abundance count at each survey time.
#' }
#' 
#' The \code{kf_survival_effects} model state provides a method for specifying the \emph{difference} in survival between the N-mixtue individuals and the 
#' known-fate individuals. Thus, the known-fate individual survival model is equal to \code{~survival_dnm + kf_survival_effects}.
#'     
#' @return A function of the vector \code{par}
#' @author Devin S. Johnson
#' @export
make_ikfdnm_likelihood = function(dnm_survival=~1, dnm_recruit=~1, dnm_det=~1, kf_survival_effects=NULL, 
                                  fixed_list=NULL, data, N_max){
  X_dnm = model.matrix(dnm_survival, data)
  if(!is.null(kf_survival_effects)){
    X_kf = cbind(X_dnm,model.matrix(kf_survival_effects, data))
  } else{
    X_kf=X_dnm
  }
  W = model.matrix(dnm_recruit, data)
  Z = model.matrix(dnm_det, data)
  np_dnm = ncol(X_dnm)
  np_kf = ncol(X_kf)
  np_rec = ncol(W)
  np_det = ncol(Z)
  new_group = as.integer(c(1,as.numeric(diff(as.numeric(data$group))!=0)))
  fixed_mat = matrix(NA, nrow=nrow(data), ncol=3)
  if(!is.null(fixed_list$survial)) fixed_mat[,1]=fixed_list$survival
  if(!is.null(fixed_list$recruit)) fixed_mat[,2]=fixed_list$recruit
  if(!is.null(fixed_list$det)) fixed_mat[,3]=fixed_list$det
  n2ll = function(par){
    beta_kf=par[1:np_kf]
    beta_dnm = par[1:np_dnm]
    rho=par[(np_kf+1):(np_kf+np_rec)]
    alpha=par[(np_kf+np_rec+1):(np_kf+np_rec+np_det)]
    omega_dnm = plogis(X_dnm%*%beta_dnm)
    omega_dnm=ifelse(is.na(fixed_mat[,1]), omega_dnm, fixed_mat[,1])
    omega_kf = plogis(X_kf%*%beta_kf)
    omega_kf=ifelse(is.na(fixed_mat[,1]), omega_kf, fixed_mat[,1])
    gamma = exp(W%*%rho)
    gamma=ifelse(is.na(fixed_mat[,2]), gamma, fixed_mat[,2])
    p=plogis(Z%*%alpha)
    p=ifelse(is.na(fixed_mat[,3]), p, fixed_mat[,3])
    out = kfdnm:::ikfdnm_hmm(
      n=data$n, 
      Y=data$Y,
      R=data$R, 
      new_group=new_group, 
      omega_dnm=omega_dnm, 
      omega_kf=omega_kf,
      gamma=gamma, 
      p=p, 
      N_max=N_max, 
      back_sample=FALSE)
    return(out[[1]]+out[[2]])
  }
  attr(n2ll,"npar") = np_kf + np_rec+np_det
  return(n2ll)
}

#' @title Create IKFDNM model MCMC sampler for abundance process
#' 
#' @description This function assimilates the data and desirted model to create a sampling function that is a function with only one argument, \code{par},
#' This likelihood function can placed within an MCMC routine for sampling the abundance process.
#' 
#' @param dnm_survival A model statement for the survival parameters of the N-mixture portion of the model
#' @param dnm_recruit A model statement for the recruitment portion of the N-mixture model
#' @param dnm_det A model statement for the the detection portion of the N-mixture model
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
#' \item \code{num_released} The number of known-fate individuals released into the population at the time of the survey
#' \item \code{num_returns} The number of known-fate individuals observed at each survey time. These are a binomial sample of the \code{num_release} of the 
#' previous row.
#' \item \code{n} The observed abundance count at each survey time.
#' }
#' 
#' The \code{kf_survival_effects} model state provides a method for specifying the \emph{difference} in survival between the N-mixtue individuals and the 
#' known-fate individuals. Thus, the known-fate individual survival model is equal to \code{~survival_dnm + kf_survival_effects}.
#'     
#' @return A function of the vector \code{par}
#' @author Devin S. Johnson
#' @export
make_ikfdnm_sampler = function(dnm_survival=~1, dnm_recruit=~1, dnm_det=~1, kf_survival_effects=NULL, 
                                  fixed_list=NULL, data, N_max){
  X_dnm = model.matrix(dnm_survival, data)
  if(!is.null(kf_survival_effects)){
    X_kf = cbind(X_dnm,model.matrix(kf_survival_effects, data))
  } else{
    X_kf=X_dnm
  }
  W = model.matrix(dnm_recruit, data)
  Z = model.matrix(dnm_det, data)
  np_dnm = ncol(X_dnm)
  np_kf = ncol(X_kf)
  np_rec = ncol(W)
  np_det = ncol(Z)
  new_group = as.integer(c(1,as.numeric(diff(as.numeric(data$group))!=0)))
  fixed_mat = matrix(NA, nrow=nrow(data), ncol=3)
  if(!is.null(fixed_list$survial)) fixed_mat[,1]=fixed_list$survival
  if(!is.null(fixed_list$recruit)) fixed_mat[,2]=fixed_list$recruit
  if(!is.null(fixed_list$det)) fixed_mat[,3]=fixed_list$det
  sampler = function(par){
    beta_kf=par[1:np_kf]
    beta_dnm = par[1:np_dnm]
    rho=par[(np_kf+1):(np_kf+np_rec)]
    alpha=par[(np_kf+np_rec+1):(np_kf+np_rec+np_det)]
    omega_dnm = plogis(X_dnm%*%beta_dnm)
    omega_dnm=ifelse(is.na(fixed_mat[,1]), omega_dnm, fixed_mat[,1])
    omega_kf = plogis(X_kf%*%beta_kf)
    omega_kf=ifelse(is.na(fixed_mat[,1]), omega_kf, fixed_mat[,1])
    gamma = exp(W%*%rho)
    gamma=ifelse(is.na(fixed_mat[,2]), gamma, fixed_mat[,2])
    p=plogis(Z%*%alpha)
    p=ifelse(is.na(fixed_mat[,3]), p, fixed_mat[,3])
    out = ikfdnm_hmm(
      n=data$n, 
      Y=data$num_returns,
      M=data$num_release,
      R=data$R, 
      new_group=new_group, 
      omega_dnm=omega_dnm, 
      omega_kf=omega_kf,
      gamma=gamma, 
      p=p, 
      N_max=N_max, 
      back_sample=TRUE)
    return(out)
  }
  attr(n2ll,"npar") = np_kf + np_rec+np_det
  return(sampler)
}


#' @title Perform MCMC sampling for integrated known-fate, dynamic N-mixture model
#' @description This function used a KFDNM likelihood function, as well as, a user defined prior distribution
#' for the parameters and uses it to perform MCMC sampling of the parameter vector. 
#' 
#' @param kfdnm_lik A likelihood function created with the \code{\link{make_ikfdnm_likelihood}} function
#' @param prior_dist A user defined function that evaluates to the log-prior distribution of the parameter vector
#' @param sample_N Logical indicating that the abundance process should be sampled.
#' @param sample_SG Logical indicating wheather the survival (S) and recruitement process should be sampled.
#' @param burn The number of iterations for burnin of the MCMC algorithm
#' @param iter The number of MCMC iterations
#' @param block Block size for adaptation of proposal variance. 
#' @param par Initial parameter value. If not provided a optimization is executed first to obtain 
#' realistic initial values. 
#' @return A list containg the MCMC sample of parameters, abundance, and dynamic abundance components.   

ikfdnm_mcmc = function(dnm_survival=~1, dnm_recruit=~1, dnm_det=~1, kf_survival_effects=NULL, 
                       fixed_list=NULL, data, N_max, ln_prior, sample_abundance=FALSE,
                       burn, iter, block, par){
  n2ll=make_ikfdnm_likelihood(
    dnm_survival=dnm_survival, dnm_recruit=dnm_recruit, dnm_det=dnm_det, 
    kf_survival_effects=kf_survival_effects, fixed_list=fixed_list, 
    data=data, N_max=N_max
    )
  NSG_sampler = make_ikfdnm_sampler(
    dnm_survival=dnm_survival, dnm_recruit=dnm_recruit, dnm_det=dnm_det, 
    kf_survival_effects=kf_survival_effects, fixed_list=fixed_list, 
    data=data, N_max=N_max
    )
  d = attr(kfdnm_lik, "npar")
  ln_post = function(par){-n2ll(par)/2 + prior_dist(par)}
  if(missing(par)){
    message("Computing starting values...")
    par_curr=optim(rep(0,d), ln_post, control=list(fnscale=-0.5))$par
  } else{par_curr=par}
  ln_post_curr = ln_post(par_curr)
  tune = 2.4^2/d
  Sigma = 0.01*diag(d)
  r_opt = 0.234
  r_hist = rep(0,burn+iter)
  par_store = matrix(NA, burn+iter, d)
  ### Update par
  message("Begininng MCMC...")
  flush.console()
  pb = txtProgressBar(min = 1, max = burn+iter, style = 3)
  for(i in 1:(burn+iter)){
    par_prop = par_curr + sqrt(tune)*crossprod(chol(Sigma),rnorm(d,0,1))
    ln_post_prop = ln_post(par_prop)
    mh = exp(ln_post_prop - ln_post_curr)
    if(runif(1,0,1)<=mh){
      par_curr = par_prop
      r_hist[i] = 1
      ln_post_curr =ln_post_prop
    }
    par_store[i,]=par_curr
    ### Adapt par update proposal
    if(i%%block==0){
      r_hat = sum(r_hist[(i-block+1):i])/block
      Sigma_hat = var(par_store[(i-block+1):i,])
      tune=exp(log(tune) + 2*((i/block)^-0.8)*(r_hat-r_opt))
      Sigma = Sigma + ((i/block)^-0.8)*(Sigma_hat-Sigma)
    }
    ### Update abundance
    if(sample_abundance){
      abund=NSG_sampler(par_curr)
    }
    
    setTxtProgressBar(pb,i)
  }
  close(pb)
  message("Processing posterior sample for output...")
  return(list(
    par_sample=par_store,
    N=NULL,
    S=NULL,
    G=NULL,
    accept=r_hist
    ))
}