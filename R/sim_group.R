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
  out = matrix(NA, num_years*num_surveys, 10)
  colnames(out)=c("year","survey","num_release","num_returns","recruit","group_surv","R","N","n","perfect")
  out[,1]=rep(1:num_years, each=num_surveys)
  out[,2]=rep(1:num_surveys, num_years)
  out[,10]=rbinom(num_years*num_surveys, 1, perfect_survey_rate)
  out[,5]=0
  out[,7]=0
  out[1,3]=num_kf
  out[1,5]=rpois(1,init_rate)
  out[1,8]=out[1,5]
  for(r in 2:(num_years*num_surveys)){
    out[r,4]=rbinom(1, out[r-1,3], survival_kf)
    out[r,6]=rbinom(1, out[r-1,8], survival_dnm)
    if(out[r,2]==1) {
      out[r,5]=rpois(1,recruit_rate)
      out[r,7] = min(3-out[r,4], out[r,6])
      out[r,3]=out[r,4] + out[r,7]
    } else {
      out[r,3]=out[r,4]
    }
    out[r,8]=out[r,5]+out[r,6]-out[r,7]
  }
  out[,9]=rbinom(n=nrow(out), size=out[,8], prob=detection)
  out[,9]=ifelse(out[,10]==1, out[,8], out[,9])
  return(as.data.frame(out))
}
