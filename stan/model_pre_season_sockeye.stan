data{
  int S; // number of periods in season
  int Y; // number of years
  int sites; // number of sites
  int netsets; // max number of net sets
  int obs; //number of observations
  int day[obs];//day of season for obervations
  int site[obs]; // site for each obs
  int set_num[obs]; // set number for each obs
  int period[obs];//day of season for obervations
  vector[obs] effort;// effort in hrs
  vector[obs] log_trs_anom;//log total runsize anomoly for year corresponding to observation
  vector[obs] zl_flow;//log total runsize anomoly for year corresponding to observation
  //vector[obs] lm_trs log(mean(total run size)) 
  int<lower=0> SockeyeC[obs];//sockeye catch
  int year[obs];// year index for obs
  int Y_obs; //years we are making estimates for (do no include any years without test fishery data)
  int ylist[Y_obs]; ///list of unique year indexes 

  vector[Y_obs-1] mu_log_forcast_anom;// log of preaseason forecast median minus lm_trs
  vector[Y_obs-1] sd_log_forcast_anom;// sd(log(preseason forecast pdf minus lm_trs))
}
parameters{
  real<lower=0> sigma_disp;
  real<lower=0> sigma_day;
  real<lower=0> sigma_set;
  real<lower=0> sigma_period;
  vector[S-1] eps_day; 
  vector[netsets-1] eps_set;  
  vector[(S*Y)-1] eps_period;
  real b4;
  real intercept;
  real b3_est;
}
transformed parameters{
  vector[S] resid_day;
  vector[S*Y] resid_period;
  vector[netsets] b2;
  vector[obs] log_lambda;
  vector[Y_obs] b3;
  resid_day[1] = 0;
  resid_period[1] = 0;
  b2[1] = 0;
  for(i in 2:S){
    resid_day[i] = resid_day[i-1] + eps_day[i-1] * sigma_day;
  }
  for(i in 2:(S*Y)){
    resid_period[i] = resid_period[i-1] + eps_period[i-1] * sigma_period;
  }
  for(i in 2:netsets){
    b2[i] =  eps_set[i-1] * sigma_set;
  }
  for(i in 1:(Y_obs-1)){
    b3[ylist[i]] = log_trs_anom[i];
  }
  b3[ylist[Y_obs]] = b3_est;
  for(i in 1:obs){
    log_lambda[i] = intercept + b2[set_num[i]] + log(effort[i]) + b3[year[i]] + resid_day[day[i]] + b4 * zl_flow[i] + resid_period[period[i]];
  }
  
}

model
{
  //priors
  sigma_day ~ std_normal();
  sigma_period ~ std_normal();
  sigma_set ~ std_normal();
  sigma_disp ~ std_normal();
  b4 ~ normal(0,5);
  intercept ~ normal(0,10);
  eps_day ~ std_normal();
  eps_set ~ std_normal();
  eps_period ~ std_normal();
  b3_est ~ normal(mu_log_forcast_anom,sd_log_forcast_anom);
  //Likelihood
  SockeyeC ~ neg_binomial_2_log(log_lambda, 1/square(sigma_disp));
}
