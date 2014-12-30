require(kfdnm)

###
### Simulation parameters
###
reps = 20
num_groups=3
num_kf = 3
num_years = 5 
num_surveys = 10
recruit_rate = 4 
init_rate = 6 
annual_survival_dnm = 0.6 
annual_survival_kf = 0.8 
perfect_survey_rate = 0.1 
detection = 0.5

set.seed(111)

par_store=matrix(NA, nrow=reps, 4)
se_store=matrix(NA, nrow=reps, 4)

###
### Begin simulation
###
pb <- txtProgressBar(min = 1, max = reps, style = 3)
for(r in 1:reps){
  data=NULL
  s_kf=annual_survival_kf^(1/num_surveys)
  s_dnm=annual_survival_dnm^(1/num_surveys)
  for(i in 1:num_groups){
    data = rbind(data, sim_group(num_kf, num_years, num_surveys, recruit_rate, init_rate, survival_dnm=s_dnm, 
                           survival_kf=s_kf, perfect_survey_rate, detection))
  }
  data$group = rep(1:num_groups, each=num_years*num_surveys)
  
  ###
  ### Create list of fixed parameters
  ###
  
  fixed_list = list(survival=rep(NA, nrow(data)), 
                    recruit=ifelse(data$year!=1 & data$survey==1, NA, 0),
                    det=ifelse(data$perfect==1, 1, NA)
  )
  
  ###
  ### Make likelihood function
  ### 
  n2ll_kfdnm = make_ikfdnm_likelihood(dnm_survival=~1, dnm_recruit=~1, dnm_det=~1, kf_survival_effects=~1, 
                                fixed_list=fixed_list, data=data, N_max=50)
  
  ###
  ### Optimize and obtain estimates and variance-covariance matrix
  ### 
  
  par_start=c(2.948662, 0.8426851, log(4), qlogis(0.5))
  n2ll_kfdnm(par_start)
  mle=optim(par_start, n2ll_kfdnm, method="BFGS", hessian=TRUE)
  par_store[r,]=mle$par
  se_store[r,] = sqrt(diag(2*solve(mle$hessian)))
  setTxtProgressBar(pb, r)
}
close(pb)
save(list=ls(), file="sim_sep_surv.RData")

