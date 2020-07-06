//=========================================================================== //
//                            Meta-regression analysis 
// "generated quantities " block to produce predicted probabilities 
// of becoming R carrier per day as a function of abx duration. 
//=========================================================================== //

// Data taken from the define_data.R file 
data {
  int<lower=0> T;                           // number of trials 
  int cases[T, 2];                          // number of outcomes reported
  int denoms[T, 2];                         // element [i,j] gives denominators for trial i arm j
  int followupdays[T];                      // number of days individuals were followed up for
  real abx_dur[T, 2];                        // abx dur is for trial i arm j
  int samples[T, 2];                        // number of surveillence samples 
}

// Parameters estimated from the model 
parameters {
  vector[15] a;                                // intercepts
  vector[15] b;                                // slopes abx duration
  real a0;                                    // mean intercept
  real b0;                                    // mean abx duration slope
  real<lower=0> sigmasq_a;                    // intercept variance 
  real<lower=0> sigmasq_b;                    // slope variance abx duration
}

transformed parameters {
  
  real sigma_a = sqrt(sigmasq_a);
  real sigma_b = sqrt(sigmasq_b);
  
  // Probability of becoming R carrier PER DAY: allow different probabilities for each arm of each trial
  real p1_1;  //trial one arm 1 
  real p1_2;  //trial one arm 2
  real p2_1;  //trial two arm 1
  real p2_2;  //trial two arm 2
  real p3_1;  //trial three arm 1
  real p3_2;  //trial three arm 2
  real p4_1;  //trial four arm 1
  real p4_2;  //trial four arm 2
  real p5_1;  
  real p5_2; 
  real p6_1;  
  real p6_2;  
  real p7_1;   
  real p7_2;  
  real p8_1;  
  real p8_2; 
  real p9_1;  
  real p9_2;  
  real p10_1;  
  real p10_2; 
  real p11_1;  
  real p11_2;  
  real p12_1;  
  real p12_2;  
  real p13_1;  
  real p13_2;  
  real p14_1;  
  real p14_2;  
  real p15_1;  
  real p15_2;  
  
  // Probabilities of outcome over FOLLOW UP PERIOD
  
  // For trials other than Wittekamp, de Smet, Hoberman,  
  // reported one outcome per participant. 
  // Hence the probabilities of outcome during the follow up period 
  // are given by 1-(1- probability of infection per day)^followupdays. 

  real pi2_1;
  real pi2_2;
  
  real pi4_1;
  real pi4_2;
  
  real pi5_1;  
  real pi5_2;
  
  real pi6_1;  
  real pi6_2;
  
  real pi7_1;  
  real pi7_2;
  
  real pi8_1;  
  real pi8_2;
  
  real pi9_1;  
  real pi9_2;
  
  real pi10_1;  
  real pi10_2;
  
  real pi11_1;  
  real pi11_2;
  
  real pi12_1;  
  real pi12_2;
  
  real pi13_1;  
  real pi13_2;
  
  real pi15_1;  
  real pi15_2;
  
  p2_1 = inv_logit(a[2] + b[2]*abx_dur[2,1]);
  p2_2 = inv_logit(a[2] + b[2]*abx_dur[2,2]);
  pi2_1 = 1-(1-p2_1)^followupdays[2];
  pi2_2 = 1-(1-p2_2)^followupdays[2];
  
  p4_1 = inv_logit(a[4] + b[4]*abx_dur[4,1]);
  p4_2 = inv_logit(a[4] + b[4]*abx_dur[4,2]);
  pi4_1 = 1-(1-p4_1)^followupdays[4];
  pi4_2 = 1-(1-p4_2)^followupdays[4];
  
  p5_1 = inv_logit(a[5] + b[5]*abx_dur[5,1]);
  p5_2 = inv_logit(a[5] + b[5]*abx_dur[5,2]);
  pi5_1 = 1-(1-p5_1)^followupdays[5];
  pi5_2 = 1-(1-p5_2)^followupdays[5];
  
  p6_1 = inv_logit(a[6] + b[6]*abx_dur[6,1]);
  p6_2 = inv_logit(a[6] + b[6]*abx_dur[6,2]);
  pi6_1 = 1-(1-p6_1)^followupdays[6];
  pi6_2 = 1-(1-p6_2)^followupdays[6];
  
  p7_1 = inv_logit(a[7] + b[7]*abx_dur[7,1]);
  p7_2 = inv_logit(a[7] + b[7]*abx_dur[7,2]);
  pi7_1 = 1-(1-p7_1)^followupdays[7];
  pi7_2 = 1-(1-p7_2)^followupdays[7];
  
  p8_1 = inv_logit(a[8] + b[8]*abx_dur[8,1]);
  p8_2 = inv_logit(a[8] + b[8]*abx_dur[8,2]);
  pi8_1 = 1-(1-p8_1)^followupdays[8];
  pi8_2 = 1-(1-p8_2)^followupdays[8];
  
  p9_1 = inv_logit(a[9] + b[9]*abx_dur[9,1]);
  p9_2 = inv_logit(a[9] + b[9]*abx_dur[9,2]);
  pi9_1 = 1-(1-p9_1)^followupdays[9];
  pi9_2 = 1-(1-p9_2)^followupdays[9];
  
  p10_1 = inv_logit(a[10] + b[10]*abx_dur[10,1]);
  p10_2 = inv_logit(a[10] + b[10]*abx_dur[10,2]);
  pi10_1 = 1-(1-p10_1)^followupdays[10];
  pi10_2 = 1-(1-p10_2)^followupdays[10];
  
  p11_1 = inv_logit(a[11] + b[11]*abx_dur[11,1]);
  p11_2 = inv_logit(a[11] + b[11]*abx_dur[11,2]);
  pi11_1 = 1-(1-p11_1)^followupdays[11];
  pi11_2 = 1-(1-p11_2)^followupdays[11];
  
  p12_1 = inv_logit(a[12] + b[12]*abx_dur[12,1]);
  p12_2 = inv_logit(a[12] + b[12]*abx_dur[12,2]);
  pi12_1 = 1-(1-p12_1)^followupdays[12];
  pi12_2 = 1-(1-p12_2)^followupdays[12];
  
  p13_1 = inv_logit(a[13] + b[13]*abx_dur[13,1]);
  p13_2 = inv_logit(a[13] + b[13]*abx_dur[13,2]);
  pi13_1 = 1-(1-p13_1)^followupdays[13];
  pi13_2 = 1-(1-p13_2)^followupdays[13];
  
  p15_1 = inv_logit(a[15] + b[15]*abx_dur[15,1]);
  p15_2 = inv_logit(a[15] + b[15]*abx_dur[15,2]);
  pi15_1 = 1-(1-p15_1)^followupdays[15];
  pi15_2 = 1-(1-p15_2)^followupdays[15];
  
  // For Wittekamp, de Smet, Hoberman,  
  // these trials reported multiple outcomes per participants. 
  // Hence the probabilities of outcome during the follow up period 
  // can be directly obtained frome the models below
  
  p1_1 = inv_logit(a[1] + b[1]*abx_dur[1,1]);
  p1_2 = inv_logit(a[1] + b[1]*abx_dur[1,2]);
  
  p3_1 = inv_logit(a[3] + b[3]*abx_dur[3,1]);
  p3_2 = inv_logit(a[3] + b[3]*abx_dur[3,2]);
  
  p14_1 = inv_logit(a[14] + b[14]*abx_dur[14,1]);
  p14_2 = inv_logit(a[14] + b[14]*abx_dur[14,2]);
  
}

model { 
  
  // For trials other than Wittekamp, de Smet, Hoberman,  
  // reported one outcome per participant 
  
  cases[2,1] ~ binomial_logit(denoms[2,1], pi2_1);
  cases[2,2] ~ binomial_logit(denoms[2,2], pi2_2);
  
  cases[4,1] ~ binomial_logit(denoms[4,1], pi4_1);
  cases[4,2] ~ binomial_logit(denoms[4,2], pi4_2);
  
  cases[5,1] ~ binomial_logit(denoms[5,1], pi5_1);
  cases[5,2] ~ binomial_logit(denoms[5,2], pi5_2);
  
  cases[6,1] ~ binomial_logit(denoms[6,1], pi6_1);
  cases[6,2] ~ binomial_logit(denoms[6,2], pi6_2);
  
  cases[7,1] ~ binomial_logit(denoms[7,1], pi7_1);
  cases[7,2] ~ binomial_logit(denoms[7,2], pi7_2);
  
  cases[8,1] ~ binomial_logit(denoms[8,1], pi8_1);
  cases[8,2] ~ binomial_logit(denoms[8,2], pi8_2);
  
  cases[9,1] ~ binomial_logit(denoms[9,1], pi9_1);
  cases[9,2] ~ binomial_logit(denoms[9,2], pi9_2);
  
  cases[10,1] ~ binomial_logit(denoms[10,1], pi10_1);
  cases[10,2] ~ binomial_logit(denoms[10,2], pi10_2);
  
  cases[11,1] ~ binomial_logit(denoms[11,1], pi11_1);
  cases[11,2] ~ binomial_logit(denoms[11,2], pi11_2);
  
  cases[12,1] ~ binomial_logit(denoms[12,1], pi12_1);
  cases[12,2] ~ binomial_logit(denoms[12,2], pi12_2);
  
  cases[13,1] ~ binomial_logit(denoms[13,1], pi13_1);
  cases[13,2] ~ binomial_logit(denoms[13,2], pi13_2);
  
  cases[15,1] ~ binomial_logit(denoms[15,1], pi15_1);
  cases[15,2] ~ binomial_logit(denoms[15,2], pi15_2);
  
  // For Wittekamp, de Smet, Hoberman,  
  // these trials reported multiple outcomes per participant. 
  
  cases[1,1] ~ binomial_logit(samples[1,1], p1_1);
  cases[1,2] ~ binomial_logit(samples[1,2], p1_2);
  
  cases[3,1] ~ binomial_logit(samples[3,1], p3_1);
  cases[3,2] ~ binomial_logit(samples[3,2], p3_2);
  
  cases[14,1] ~ binomial_logit(samples[14,1], p14_1);
  cases[14,2] ~ binomial_logit(samples[14,2], p14_2);
  
  a0 ~ normal(0, 10);
  b0 ~ normal(0, 10);
  
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
  
  a_new = normal_rng(a0, sigmasq_a );
  b_new = normal_rng(b0, sigmasq_b );
  
  for (i in 1:14){
    // calculate mean predicted probability of becoming a R carrier as a function of abx days
    predicted_mean[i] =  inv_logit(a0 + b0*i);
  }  
  
}  
