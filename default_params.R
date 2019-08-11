############## Shared parameters ##################
#fixed parameters 
n.bed=30                            # n.bed= number of beds in the ward
mean.max.los=7                      # mean.max.los= mean of max length of stay (exponential distribution)
n.day=365

############# Simple model #################
short_dur = 5
long_dur = 14
meanDur=10

p.infect = 0.2                   # p=probability of receiving antibiotic
p.r.day1 =0.1
prob_StartBact_R = 0.2           # Probability of being colonized with resistant strain on admission
prop_S_nonR = 0.8                # Proportion of large S within non-resistant states (S+s)
bif = 0.9                        # bacterial interference factor 
pi_ssr = 0.01                   # pi_s= probability of R transmitting to ss
repop.s1 = 0.01
mu_r = 0.02                     # mu_r= rate of clearance of R to become S 
abx.s = 0.5                      # probability of clearing S to become ss under antibiotic treatment 
abx.r = 0.5
cum.r.1=1000

############## Binary model ###################
prop_S_nonR = 0.8                 # Proportion of large S within non-resistant states (S+s)
prop_Sr_inR = 0.7                  
prop_sr_inR = 0.1                     
mu1 = 0.015                       # mu1= probability of clearance of Sr to become S
mu2 = 0.015                       # mu2= probability of clearance of sr to become s 
mu_r=0.015 
repop.r1= 0.01                    # probability of repopulation of Sr to become sR 
repop.r2= 0.01
repop.s1= 0.005                   # probability of repopulation of s to become S 
repop.s2= 0.005                   # probability of repopulation of sr to become SR 

############# Frequency model #################
K = 11                             # gut holding capacity
total_prop = 0.02                   # mean of total starting amount of gut bacteria on log scale
#t_sd = 1.1218                       # sd of total starting amount of gut bacteria on log scale
r_prop = 0.3                        # mean of starting amount of resistant gut bacteria on log scale
#r_sd = 1.8921                       # sd of starting amount of resistant gut bacteria on log scale
r_growth = 0.4                        # r_growth = growth constant for logistic growth
s_growth = 0.3
r_thres = 5                        # r_thres = R threshold level for tranmissibility
r_trans = 6                        # r_trans = amount transmitted
abx.s = 10                          # abxr_killr = amount of r killed by broad spectrum abx r
abx.r = 10                          # abxs_kills = amount of s killed by narrow spectrum abx s


if (model=='simple' | model=='binary'){
    abx.s = 0.5                      # probability of clearing S to become ss under antibiotic treatment 
    abx.r = 0.5
} else {
    abx.s = 10                       # abxr_killr = amount of r killed by broad spectrum abx r
    abx.r = 10                       # abxs_kills = amount of s killed by narrow spectrum abx s
}
    
