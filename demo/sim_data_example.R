#require(kfdnm)

###
### Create data to mimic wolf pack example
###

set.seed(111)
data=NULL
for(i in 1:20){
  data = rbind(data,
              sim_group(
                num_KF = 3, 
                num_years = 5, 
                num_surveys = 10, 
                recruit_rate = 4, 
                init_rate = 6, 
                annual_survival = 0.6, 
                annual_survival_KF = 0.8, 
                perfect_survey_rate = 0.1, 
                detection = 0.5
              )
  )
}
data$group = rep(1:20, each=5*10)

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
llkf = make_ikfdnm_likelihood(dnm_survival=~1, dnm_recruit=~1, dnm_det=~1, kf_survival_effects=~1, 
                                  fixed_list=fixed_list, data=data, N_max=50)

###
### Optimize and obtain estimates and variance-covariance matrix
### 

par_start=c(qlogis(0.6^0.1), qlogis(0.8^0.1)-qlogis(0.6^0.1), log(4), qlogis(0.5))

mle=optim(par_start, llkf, method="BFGS", control=list(REPORT=1, trace=1), hessian=TRUE)
par=mle$par
V = 2*solve(mle$hessian)


###
### Estimates and 95% CI
###

# omega_dnm
cat("N-mixture annual survival: \n")
cat(plogis(par[1])^10, "(", plogis(par[1]-2*sqrt(V[1,1]))^10, ",",plogis(par[1]+2*sqrt(V[1,1]))^10, ")\n")

# omega_kf
se = sqrt(sum(V[1:2,1:2]))
cat("Known-fate annual survival: \n")
cat(plogis(par[1]+par[2])^10, "(", plogis(par[1]+par[2]-2*se)^10, ",",plogis(par[1]+par[2]+2*se)^10, ")\n")

# gamma
cat("Recruitment rate:")
cat(exp(par[3]), "(", exp(par[3]-2*sqrt(V[3,3])), ",",exp(par[3]+2*sqrt(V[3,3])), ")\n")

# p
cat("Detection prob.:\n")
cat(plogis(par[4]), "(", plogis(par[4]-2*sqrt(V[4,4])), ",",plogis(par[4]+2*sqrt(V[4,4])), ")\n")

