###########################1. Input parameters##########################
#fixed parameters 
n.bed<-20                             # n.bed= number of beds in the ward
n.days<- 100                          # n.days= number of days we want to observe
mean.max.los<-20                      # mean.max.los= mean of max length of stay (normal distribution)

#variable parameters 
###epidemiological 
p.s<-0.5                              # p=probability of receiving antibiotic for sensitive organisms
p.r<-0.1                              # p= daily probability of contracting HAI and receiving antibiotic for resistant organisms 
#prob_StartBact_bi<-c(0.5,0.2, 0.1, 0.05) # prob_StartBact= vector of probability of carrying c(S,Sr, sR, sr)
#                                     # possible states: S- carry sensitive organism only 
#                                                        Sr- carry largely sensitive organisms and small population of resistant organisms
#                                                        sR- carry largely resistant organisms and small population of sensitive organisms
#                                                        sr- carry small population of sensitive organisms and resistant organisms
#                                                        s - carry small population of sensitive organisms 
prob_StartBact_R <- 0.35
prop_S_nonR <- 0.5                # Proportion of large S within non-resistant states (S+s)
prop_Sr_inR <- 0.55                      # Proportion of large S within non-resistant states (S+s)
prop_sr_inR <- 0.15                      # Proportion of large S within non-resistant states (S+s)

###biological 
pi_r1 <- 0.003                        # pi_r1= probability of R transmitting to S to become Sr
bif_bi<- 2                               # bacterial interference factor 
pi_r2 <- pi_r1 * bif_bi                  # pi_r2= probability of R transmitting to s to become sr 
#                                       (pi_r1 < pi_r2 if being colonised with S protects colonisation by R)

mu1 <- 0                              # mu1= probability of clearance of Sr to become S
mu2 <- 0                              # mu2= probability of clearance of sr to become s 

abx.s<-0.2                            # probability of clearing S to become s under antibiotic treatment (daily)
abx.r<-0.3                            # probability of clearing R to become r under antibiotic treatment (daily)

repop.s1<- 0                          # probability of repopulation of s to become S 
repop.s2<- 0                          # probability of repopulation of sr to become SR 
repop.s3<- 0.01                       # probability of repopulation of sR to become sr
repop.r1<- 0                          # probability of repopulation of Sr to become sR 
repop.r2<- 0                          # probability of repopulation of sr to become sR 
#bif1<- 0.2                             # bacterial interference factor - how much antibiotics selects for R
# repop.r3 <- repop.r2*bif1              # probability of repopulation of sr to become sR
#                                        ( repop.r3 < repop.r2 if antibiotics increases selection for R )

### in-host gut (freq model only)
bact_slots <- 1000                      # environmental carrying capacity
bact_start <- 500                       # bacteria level (number) for starting 
                                        # where big letter = 1*bact_start, singly small = 0.5*bact_start
                                        # small next to big = 0.05*bact_start
R_thres <- 100                          # R threshold level for tranmissibility
abxr_killr <- 500                       # amount of r killed by broad spectrum abx r
abxr_kills <- 500                       # amount of s killed by broad spectrum abx r
abxs_kills <- 500                       # amount of s killed by narrow spectrum abx s
r_trans <- 100                          # amount transmitted
r_growth <- 2                           # growth constant for logistic growth

mean_dur <- 4                           # antibiotic duration (days)

# Simple model has prob_Start has less slots than others

minDur<- 1                            # minDur=the min duration of antibiotic

#variable parameters 
###epidemiological 
p <- 0.2                                # p=probability of receiving antibiotic
#prob_StartBact<-c(0.4,0.58)           # prob_StartBact= vector of probability of carrying sensitive organism, resistant organism

###biological 
pi_s <- 0.003                         # pi_s= probability of S transmitting to N 
pi_r <- 0.003                         # pi_r= probability of R transmitting to N 
bif<- 0.002                               # bacterial interference factor 
pi_sr <- pi_r * bif                   # pi_sr= probability of R transmitting to S (a proportion of pi_r if being colonised with S protects colonisation by R)
mu_r <- 0                             # mu_r= rate of clearance of R to become S 
abx.clear<-0.2                        # probability of clearing S to become N under antibiotic treatment 

###gut parameters
K <- 1000
t_mean <- 4.0826
t_sd <- 1.1218
r_mean <- 1.7031
r_sd <- 1.8921
