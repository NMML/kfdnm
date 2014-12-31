require(kfdnm)

###
### Create data to mimic wolf pack example
###
num_groups=5
dnm_data = NULL
kf_data = NULL
for(i in 1:num_groups){
  data_list = sim_data(
    num_kf = 3, 
    num_years = 5, 
    num_surveys = 10, 
    recruit_rate = 4, 
    init_rate = 6, 
    survival_dnm = 0.95, 
    survival_kf = 0.95, 
    perfect_survey_rate = 0.1, 
    detection = 0.5
  )
  data_list$kf_data$id = paste(i, data_list$kf_data$id, sep="-")
  data_list$dnm_data$group=i
  dnm_data = rbind(dnm_data, data_list$dnm_data) 
  kf_data = rbind(kf_data, data_list$kf_data)
}
kf_data$id = factor(kf_data$id)
dnm_data$group = factor(dnm_data$group)


###
### Create list of fixed parameters
###

fixed_list = list(
  recruit=ifelse(dnm_data$year!=1 & dnm_data$survey==1, NA, 0),
  detection=ifelse(dnm_data$perfect==1, 1, NA)
)

###
### Make likelihood function
### 
n2ll_kfdnm = make_ikfdnm_likelihood(survival=~1, recruit=~1, detection=~1, kf_survival_effects=NULL, 
                                    fixed_list=fixed_list, kf_data=kf_data, dnm_data=dnm_data, N_max=50)

###
### Optimize and obtain estimates and variance-covariance matrix
### 

par_start=c(qlogis(0.95), log(4), qlogis(0.5))
#par_start=rep(0,3)
n2ll_kfdnm(par_start)

mle=optim(par_start, n2ll_kfdnm, method="BFGS", control=list(REPORT=1, trace=1), hessian=TRUE)
par=mle$par
se = sqrt(diag(2*solve(mle$hessian)))


###
### Estimates and 95% CI
###

# omega
cat("Survival: \n")
cat(plogis(par[1]), "(", plogis(par[1]-2*se[1]), ",",plogis(par[1]+2*se[1]), ")\n")

# gamma
cat("Recruitment rate:")
cat(exp(par[2]), "(", exp(par[2]-2*se[2]), ",",exp(par[2]+2*se[2]), ")\n")

# p
cat("Detection prob.:\n")
cat(plogis(par[3]), "(", plogis(par[3]-2*se[3]), ",",plogis(par[3]+2*se[3]), ")\n")

