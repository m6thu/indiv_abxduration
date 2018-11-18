############## Shared parameters ##################
#fixed parameters 
n.bed<-20                             # n.bed= number of beds in the ward
n.days<- 100                          # n.days= number of days we want to observe
mean.max.los<-20                      # mean.max.los= mean of max length of stay (normal distribution)
meanDur <- 7                          # meanDur = mean duration of antibiotics

############# Simple model #################
p <- 0.2                                # p=probability of receiving antibiotic
prob_StartBact_R <- 0.35             # Probability of being colonized with resistant strain on admission
prop_S_nonR <- 0.5                # Proportion of large S within non-resistant states (S+s)

pi_r <- 0.003                         # pi_r= probability of R transmitting to N 
bif<- 0.002                               # bacterial interference factor 
mu_r <- 0                             # mu_r= rate of clearance of R to become S 
abx.clear<-0.2                        # probability of clearing S to become N under antibiotic treatment 

############## Binary model ###################
p.s<-0.5                              # p=probability of receiving antibiotic for sensitive organisms
p.r.day1<-0.1                         # p= probability of receiving antibiotic for resistant organisms on day 1 admission 
p.r.dayafter<-0.1                     # p= daily probability of contracting HAI and receiving antibiotic for resistant organisms

prob_StartBact_R <- 0.35             # Probability of being colonized with resistant strain on admission
prop_S_nonR <- 0.5                # Proportion of large S within non-resistant states (S+s)
prop_Sr_inR <- 0.55                      # Proportion of large S within non-resistant states (S+s)
prop_sr_inR <- 0.15                      # Proportion of large S within non-resistant states (S+s)

pi_r1 <- 0.003                        # pi_r1= probability of R transmitting to S to become Sr
bif<- 2                               # bacterial interference factor 
mu1 <- 0                              # mu1= probability of clearance of Sr to become S
mu2 <- 0                              # mu2= probability of clearance of sr to become s 
repop.r1<- 0                          # probability of repopulation of Sr to become sR 
repop.r2<- 0                          # probability of repopulation of sr to become sR 
repop.s1<- 0                          # probability of repopulation of s to become S 
repop.s2<- 0                          # probability of repopulation of sr to become SR 
repop.s3<- 0.01                       # probability of repopulation of sR to become sr

abx.s<-0.2                            # probability of clearing S to become s under antibiotic treatment (daily)
abx.r<-0.3                            # probability of clearing R to become r under antibiotic treatment (daily)

############# Frequency model #################

# Same as binary model with gut model parameters as below
K <- 1000
t_mean <- 4.0826
t_sd <- 1.1218
r_mean <- 1.7031
r_sd <- 1.8921

r_thres <- 100                          # R threshold level for tranmissibility
r_trans <- 100                          # amount transmitted
r_growth <- 2                           # growth constant for logistic growth
abxr_killr <- 500                       # amount of r killed by broad spectrum abx r
abxr_kills <- 500                       # amount of s killed by broad spectrum abx r
abxs_kills <- 500                       # amount of s killed by narrow spectrum abx s

