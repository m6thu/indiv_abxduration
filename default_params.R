###########################1. Input parameters##########################
#fixed parameters 
n.bed<-20                             # n.bed= number of beds in the ward
n.days<- 100                          # n.days= number of days we want to observe
mean.max.los<-20                      # mean.max.los= mean of max length of stay (normal distribution)

#variable parameters 
###epidemiological 
p.s<-0.5                              # p=probability of receiving antibiotic for sensitive organisms
p.r<-0.1                              # p= daily probability of contracting HAI and receiving antibiotic for resistant organisms 
prob_StartBact<-c(0.5,0.2, 0.1, 0.05) # prob_StartBact= vector of probability of carrying c(S,Sr, sR, sr)
#                                     # possible states: S- carry sensitive organism only 
#                                                        Sr- carry largely sensitive organisms and small population of resistant organisms
#                                                        sR- carry largely resistant organisms and small population of sensitive organisms
#                                                        sr- carry small population of sensitive organisms and resistant organisms
#                                                        s - carry small population of sensitive organisms 

###biological 
pi_r1 <- 0.003                        # pi_r1= probability of R transmitting to S to become Sr
bif<- 2                               # bacterial interference factor 
pi_r2 <- pi_r1 * bif                  # pi_r2= probability of R transmitting to s to become sr 
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
bif1<- 0.2                             # bacterial interference factor - how much antibiotics selects for R
repop.r3 <- repop.r2*bif1              # probability of repopulation of sr to become sR
#                                        ( repop.r3 < repop.r2 if antibiotics increases selection for R )

mean_dur <- 4                           # antibiotic duration
