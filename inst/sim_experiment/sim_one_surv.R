require(kfdnm)

###
### Simulation parameters
###
reps = 200
num_groups=20
num_kf = 3
num_years = 5 
num_surveys = 10
recruit_rate = 4 
init_rate = 6 
annual_survival_dnm = 0.8 
annual_survival_kf = 0.8 
perfect_survey_rate = 0.1 
detection = 0.5

set.seed(111)

par_store=matrix(NA, nrow=reps, 3)
se_store=matrix(NA, nrow=reps, 3)

###
### Begin simulation
###
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
  llkf = make_ikfdnm_likelihood(dnm_survival=~1, dnm_recruit=~1, dnm_det=~1, kf_survival_effects=NULL, 
                                fixed_list=fixed_list, data=data, N_max=50)
  
  ###
  ### Optimize and obtain estimates and variance-covariance matrix
  ### 
  
  par_start=c(qlogis(annual_survival_dnm^0.1), log(4), qlogis(0.5))
  
  mle=optim(par_start, llkf, method="BFGS", hessian=TRUE)
  par_store[r,]=mle$par
  se_store[r,] = sqrt(diag(2*solve(mle$hessian)))
  cat(r,"\n")
}
save(list=ls(), file="sim_one_surv.RData")

