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
      N_max=50, 
      back_sample=FALSE)[1:2]
    return(out[[1]]+out[[2]])
  }
  return(n2ll)
}

# ikfdnm_mcmc = function(data, burn, iter){
#   
# }