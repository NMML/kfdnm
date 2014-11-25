require(kfdnm)

###
### Create data to mimic wolf pack example
###

set.seed(111)
data=NULL
num_groups=3
for(i in 1:num_groups){
  data = rbind(data,
              sim_group(
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
  )
}
data$group = rep(1:num_groups, each=5*10)

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
kfdnm_lik = make_ikfdnm_likelihood(
  dnm_survival=~1, 
  dnm_recruit=~1, 
  dnm_det=~1, 
  kf_survival_effects=NULL, 
  fixed_list=fixed_list, 
  data=data, 
  N_max=50
  )

###
### Make prior function
###
ln_laplace=function(par){-sum(abs(par))}

###
### Optimize and obtain estimates and variance-covariance matrix
### 

par_start=c(qlogis(0.9), log(3), qlogis(0.7))

post=ikfdnm_mcmc(
  kfdnm_lik, 
  prior_dist=ln_laplace, 
  sample_N=FALSE, 
  sampleSG=FALSE, 
  burn=10, 
  iter=100, 
  block=10, 
  par=par_start)


###
### Estimates and 95% CI
###

# omega
cat("Annual survival: \n")
cat(plogis(par[1])^10, "(", plogis(par[1]-2*se[1])^10, ",",plogis(par[1]+2*se[1])^10, ")\n")

# gamma
cat("Recruitment rate:")
cat(exp(par[2]), "(", exp(par[2]-2*se[2]), ",",exp(par[2]+2*se[2]), ")\n")

# p
cat("Detection prob.:\n")
cat(plogis(par[3]), "(", plogis(par[3]-2*se[3]), ",",plogis(par[3]+2*se[3]), ")\n")

