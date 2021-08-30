############## Shared parameters ##################
#fixed parameters 
n.bed<-20                             # n.bed= number of beds in the ward
mean.max.los<-7                       # mean.max.los= mean of max length of stay (exponential distribution)
sdDur <- 1                            # sdDur = sd of duration of antibiotics
n.day<-300

############# Simple model #################
short_dur <- 5
long_dur <- 20

p <- 0.4                            # p=probability of receiving antibiotic
prob_StartBact_R <- 0.4             # Probability of being colonized with resistant strain on admission
prop_S_nonR <- 0.5                  # Proportion of large S within non-resistant states (S+s)
bif <- 0.3                          # bacterial interference factor 
pi_ssr <- 0.0001                    # pi_s= probability of R transmitting to ss
repop.s1 <- 0
mu_r <- 0                           # mu_r= rate of clearance of R to become S 

############## Binary model ###################
short_dur.s <- 5
long_dur.s <- 20
short_dur.r <- 5
long_dur.r <- 20

p.s<-0.3                              # p=probability of receiving antibiotic for sensitive organisms
p.r.day1<-0.1                         # p= probability of receiving antibiotic for resistant organisms on day 1 admission 
p.r.dayafter<-0.0005                     # p= daily probability of contracting HAI and receiving antibiotic for resistant organisms
prob_StartBact_R <- 0.4             # Probability of being colonized with resistant strain on admission
prop_S_nonR <- 0.5                   # Proportion of large S within non-resistant states (S+s)
prop_Sr_inR <- 0.55                  
prop_sr_inR <- 0.15                     
pi_r1 <- 0.0001                        # pi_r1= probability of R transmitting to S to become Sr
bif<- 0.5                               # bacterial interference factor 
mu1 <- 0                              # mu1= probability of clearance of Sr to become S
mu2 <- 0                              # mu2= probability of clearance of sr to become s 
abx.s<-0.5                            # probability of clearing S to become s under antibiotic treatment
abx.r<-0.5                            # probability of clearing R to become r under antibiotic treatment
repop.r<- 0                          # probability of repopulation of Sr to become sR 
repop.s1<- 0                          # probability of repopulation of s to become S 
repop.s2<- 0                          # probability of repopulation of sr to become SR 
depop.r<- 0.001                       # probability of repopulation of sR to become sr

############# Frequency model #################
pi_r <- 0.0001                             # p.t = daily probability of transmitting resistant E coli
K <- 15                            # gut holding capacity
t_mean <- 4.0826                    # mean of total starting amount of gut bacteria on log scale
t_sd <- 1.1218                      # sd of total starting amount of gut bacteria on log scale
r_mean <- 1.7031                    # mean of starting amount of resistant gut bacteria on log scale
r_sd <- 1.8921                       # sd of starting amount of resistant gut bacteria on log scale
r_growth <- 2                           # r_growth = growth constant for logistic growth
r_thres <- 10                          # r_thres = R threshold level for tranmissibility
r_trans <- 10                          # r_trans = amount transmitted
abxr_killr <- 10                       # abxr_killr = amount of r killed by broad spectrum abx r
abxr_kills <- 10                       # abxr_kills = amount of s killed by broad spectrum abx r
abxs_kills <- 10                       # abxs_kills = amount of s killed by narrow spectrum abx s

