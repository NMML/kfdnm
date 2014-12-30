#' @title Simulate group data for analysis with known-fate dynamic N-mixture models
#' 
#' @description This function generates simulated data for a single group in the correct form for use in the kfdnm package. In addition to the required
#' elements for estimation, additional latent information is also provided. The true abundance, recruits, and survivors are retained for user reference. 
#' 
#' @param num_kf Target numer of known-fate individuals in the group
#' @param num_years Number of years for the study
#' @param num_surveys Number of surveys within each year
#' @param recruit_rate Average recruitment rate
#' @param init_rate Initial average group size in year 1 survey 1.
#' @param survival_dnm Survival rate for non-kf individuals.
#' @param survival_kf Survival rate for KF individuals.
#' @param perfect_survey_rate The probability that all non-kf individuals are observed in a survey.
#' @param detection The detection probability for abundance surveys.
#' 
#' @export
#' 
sim_group = function(num_kf, num_years, num_surveys, recruit_rate, init_rate, 
                    survival_dnm, survival_kf, perfect_survey_rate, detection){
  out = data.frame(year=rep(1:num_years, each=num_surveys))
  out$survey=rep(1:num_surveys, num_years)
  out$perfect_survey=rbinom(num_years*num_surveys, 1, perfect_survey_rate)
  out$Y = 0
  out$R = 0
  out$R[1]=num_kf
  N = rep(NA,num_years*num_surveys) 
  N[1] = rpois(1,init_rate)+num_kf
  S = rep(NA,num_years*num_surveys)
  G = rep(0,num_years*num_surveys)
  for(r in 2:(num_years*num_surveys)){
    out$Y[r]=rbinom(1, out$Y[r-1]+out$R[r-1], survival_kf)
    S[r]=rbinom(1, N[r-1]-out$R[r-1], survival_dnm)
    if(out$survey[r]==1) {
      G[r]=rpois(1,recruit_rate)
    }
    N[r] = S[r]+ G[r]
    if(out$survey[r]==1){
      out$R[r]= min(num_kf-out$Y[r], N[r])
    }
  }
  det=rep(detection, num_years*num_surveys)
  det=ifelse(out$perfect_survey==1, 1, det)
  out$S=S
  out$G=G
  out$N=N
  out$n = rbinom(num_years*num_surveys, N-out$R, det)
  return(out)
}
