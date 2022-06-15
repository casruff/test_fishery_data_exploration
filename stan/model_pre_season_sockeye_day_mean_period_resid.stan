data{
  int S; // number of periods in season
  int Y; // number of years
  //int sites; // number of sites
  int netsets; // max number of net sets
  int obs; //number of observations
  int day[obs];//day of season for obervations
  //int site[obs]; // site for each obs
  int set_num[obs]; // set number for each obs
  //int period[obs];//day of season for obervations
  vector[obs] effort;// effort in hrs
  vector[obs] zl_flow;//log total runsize anomoly for year corresponding to observation
  //vector[obs] lm_trs log(mean(total run size)) 
  int<lower=0> SockeyeC[obs];//sockeye catch
  int year[obs];// year index for obs
  int Y_obs; //years we are making estimates for (do no include any years without test fishery data)
  int ylist[Y_obs]; ///list of unique year indexes 
  vector[Y_obs-1] log_trs_anom;//log total runsize anomoly for year corresponding to observation
  real mu_log_forcast_anom;// log of preaseason forecast median minus lm_trs
  real<lower=0> sd_log_forcast_anom;// sd(log(preseason forecast pdf minus lm_trs))
}
parameters{
  real<lower=-1,upper=1> phi1;
  real<lower=-1,upper=1> phi2;
  real<lower=0> sigma_disp;
  real<lower=0> sigma_set;
  real<lower=0> sigma_period;
  real<lower=0> sigma_init;
  real<lower=0> sigma_day;
  vector[netsets-1] eps_set;  
  vector[(S*Y)] eps_period;
  vector[S-1] eps_day;
  real b4;
  real intercept;
  real b3_anom;
}
transformed parameters{
  vector[S*Y] resid_period;
  vector[S] resid_day;
  simplex[S]day_pct[Y_obs];
  vector[netsets] b2;
  vector[obs] log_lambda;
  vector[Y_obs] b3;
  resid_period[1] = eps_period[1] * sigma_init*3;
  resid_day[1] = 0;
  b2[1] = 0;
  
  for(i in 2:(S)){
    resid_day[i] = resid_day[i-1] + eps_day[i-1] * sigma_day;
  }
  
  for(i in 2:(S+1)){
    resid_period[i] = resid_period[i-1] + eps_period[i] * sigma_init;
  }
  for(i in (S+2):(S*Y)){
    resid_period[i] = phi1 * resid_period[i-1] + phi2 * (resid_period[i-S] - resid_period[i-S-1]) + eps_period[i] * sigma_period;
  }
  for(i in 1:Y_obs){
    day_pct[i][1:S] = softmax(resid_day[1:S] + resid_period[((i-1)*S+1):(i*S)]);
  }

  for(i in 2:netsets){
    b2[i] =  eps_set[i-1] * sigma_set;
  }
  for(i in 1:(Y_obs-1)){
    b3[ylist[i]] = log_trs_anom[ylist[i]];
  }
  b3[ylist[Y_obs]] = b3_anom * sd_log_forcast_anom + mu_log_forcast_anom;
  for(i in 1:obs){
    log_lambda[i] = intercept + b2[set_num[i]] + log(effort[i]) + b3[year[i]] + b4 * zl_flow[i] + log(day_pct[year[i]][day[i]]); 
  }
  
}

model
{
  //priors
  phi1 ~ std_normal();
  phi2 ~ std_normal();
  sigma_init ~ std_normal(); 
  sigma_day ~ std_normal();
  sigma_period ~ std_normal();
  sigma_set ~ std_normal();
  sigma_disp ~ std_normal();
  b4 ~ normal(0,5);
  intercept ~ normal(0,10);
  eps_set ~ std_normal();
  eps_period ~ std_normal();
  eps_day ~ std_normal();
  b3_anom ~ std_normal();
  //Likelihood
  SockeyeC ~ neg_binomial_2_log(log_lambda, 1/square(sigma_disp));
}
