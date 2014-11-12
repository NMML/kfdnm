#' @title Create IKFDNM model likehood function
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