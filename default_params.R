############## Shared parameters ##################
#fixed parameters 
n.bed<-20                             # n.bed= number of beds in the ward
n.day<- 100                         # n.days= number of days we want to observe
mean.max.los<-20                      # mean.max.los= mean of max length of stay (normal distribution)
meanDur <- 7                          # meanDur = mean duration of antibiotics
sdDur <- 5                              # sdDur = sd of duration of antibiotics

timestep <- 1                          # timestep = (int) number of timesteps per day e.g. 24 is hourly timestep
##### All probabilities are per timestep as defined here

############# Simple model #################
p <- 0.2                                # p=probability of receiving antibiotic
prob_StartBact_R <- 0.35             # Probability of being colonized with resistant strain on admission
prop_S_nonR <- 0.5                # Proportion of large S within non-resistant states (S+s)

bif <- 0.002 # bacterial interference factor 
pi_ssr <- 0.3 # pi_s= probability of R transmitting to ss
repop.s1 <- 0
mu_r <- 0  # mu_r= rate of clearance of R to become S 
abx.clear <- 0.1  # probability of clearing S to become N under antibiotic treatment 

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

abx.s<-0.2                            # probability of clearing S to become s under antibiotic treatment
abx.r<-0.3                            # probability of clearing R to become r under antibiotic treatment

############# Frequency model #################

pi_r <- 0.06                             # p.t = daily probability of transmitting resistant E coli

###### All parameters below are on log scale
K <- 100                            # gut holding capacity
t_mean <- 4.0826
t_sd <- 1.1218
r_mean <- 1.7031
r_sd <- 1.8921

r_growth <- 2                           # r_growth = growth constant for logistic growth
r_thres <- 10                          # r_thres = R threshold level for tranmissibility
r_trans <- 10                          # r_trans = amount transmitted
abxr_killr <- 5                       # abxr_killr = amount of r killed by broad spectrum abx r
abxr_kills <- 5                       # abxr_kills = amount of s killed by broad spectrum abx r
abxs_kills <- 5                       # abxs_kills = amount of s killed by narrow spectrum abx s

