//=========================================================================== //
//                            Meta-regression analysis 
// "generated quantities" block to produce predicted probabilities 
// of becoming R carrier per day as a function of abx duration. 
//=========================================================================== //


// Data taken from the define_data.R file 
data {
  int<lower=0> T;                           // number of trials 
  int cases[T, 2];                          // number of outcomes reported
  int denoms[T, 2];                         // element [i,j] gives denominators for trial i arm j
  int followupdays[T];                      // number of days individuals were followed up for
  real abx_dur[T, 2];                       // abx dur is for trial i arm j
  int samples[T, 2];                        // number of surveillence samples 
}


// Constant unmodeled parameters
transformed data {
  int N = 15;                                               // number of trials
  int M = 2;                                                // number of arms per trial
  int multi_outcomes[3] = {1,3,14};                         // indices of multiple outcome trials (i.e. Wittekamp, de Smet, Hoberman)
  int single_outcomes[12] = {2,4,5,6,7,8,9,10,11,12,13,15}; // indices of single   outcome trials (i.e. the others)
}


// Parameters estimated from the model 
parameters {
  vector[15] a;                             // intercepts
  vector[15] b;                             // slopes abx duration
  real a0;                                  // mean intercept
  real b0;                                  // mean abx duration slope
  real<lower=0> sigmasq_a;                  // intercept variance 
  real<lower=0> sigmasq_b;                  // slope variance abx duration
}


transformed parameters {
  real sigma_a = sqrt(sigmasq_a);
  real sigma_b = sqrt(sigmasq_b);
  
  real p[N,M];         // Probability of becoming R carrier PER DAY      (p[i,j]        = probability for trial i, arm j)
  real p_follow[N,M];  // Probabilities of outcome over FOLLOW UP PERIOD (p_follow[i,j] = probability for trial i, arm j)
  
  for (i in 1:N) {
    for (j in 1:M) {
      p[i,j] = a[i] + b[i]*abx_dur[i,j];
      p_follow[i,j] = 1-(1-p[i,j])^followupdays[i];
    }
  }
}


model {
  // Wittekamp, de Smet, Hoberman reported *multiple* outcomes per participants. 
  // Hence the probabilities of outcome during the follow up period can be directly obtained frome the model (i.e. `p`)
  for (i in multi_outcomes) {
    for (j in 1:M) {
      cases[i,j] ~ binomial_logit(samples[i,j], p[i,j]);
    }
  }
  
  // The other trials reported *one* outcome per participant. 
  // Hence the probabilities of outcome during the follow up period are given by 
  // `1 - (1 - probability of infection per day)^followupdays` (i.e. `pii`)
  for (i in single_outcomes) {
    for (j in 1:M) {
      cases[i,j] ~ binomial_logit(denoms[i,j], p_follow[i,j]);
    }
  }
  
  a0 ~ normal(1, 10);
  b0 ~ normal(1, 10);
  
  a ~ normal(a0, sigmasq_a);
  b ~ normal(b0, sigmasq_b);
  
  sigmasq_a ~ normal(0, 1); 
  sigmasq_b ~ normal(0, 1); 
}


generated quantities {
  real OR_abx_dur = exp(b0);
  
  real a_new;
  real b_new;
  
  real<lower=0, upper=1> predicted_mean[14];
  
  a_new = normal_rng(a0, sigmasq_a);
  b_new = normal_rng(b0, sigmasq_b);
  
  for (i in 1:14){
    // calculate mean predicted probability of becoming a R carrier as a function of abx days
    predicted_mean[i] =  inv_logit(a0 + b0*i);
  }  
}  
