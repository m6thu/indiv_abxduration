############## Shared parameters ##################
#fixed parameters 
n.bed=30                            # n.bed= number of beds in the ward
mean.max.los=7                      # mean.max.los= mean of max length of stay (exponential distribution)
n.day=440

############# Simple model #################
short_dur = 5
long_dur = 21

p.infect = 0.4                   # p=probability of receiving antibiotic
p.r.day1 =0.2
prob_StartBact_R = 0.4           # Probability of being colonized with resistant strain on admission
prop_S_nonR = 0.8                # Proportion of large S within non-resistant states (S+s)
bif = 0.3                        # bacterial interference factor 
pi_ssr = 0.005                   # pi_s= probability of R transmitting to ss
repop.s1 = 0.001
mu_r = 0.001                     # mu_r= rate of clearance of R to become S 
abx.s = 0.5                      # probability of clearing S to become ss under antibiotic treatment 
abx.r = 0.5
cum.r.1=1000

############## Binary model ###################
prob_StartBact_R = 0.4               # Probability of being colonized with resistant strain on admission
prop_S_nonR = 0.8                    # Proportion of large S within non-resistant states (S+s)
prop_Sr_inR = 0.8                  
prop_sr_inR = 0.1                     
mu1 = 0.001                          # mu1= probability of clearance of Sr to become S
mu2 = 0.001                          # mu2= probability of clearance of sr to become s 
mu_r=0.001
repop.r1= 0.001                      # probability of repopulation of Sr to become sR 
repop.r2= 0.001
repop.s1= 0.001                      # probability of repopulation of s to become S 
repop.s2= 0.001                      # probability of repopulation of sr to become SR 

############# Frequency model #################
K = 15                            # gut holding capacity
total_prop = 4.0826/15      # mean of total starting amount of gut bacteria on log scale
t_sd = 1.1218                      # sd of total starting amount of gut bacteria on log scale
r_prop = 1.7031/4.0826          # mean of starting amount of resistant gut bacteria on log scale
r_sd = 1.8921                       # sd of starting amount of resistant gut bacteria on log scale
r_growth = 2                           # r_growth = growth constant for logistic growth
s_growth = 2
r_thres = 10                          # r_thres = R threshold level for tranmissibility
r_trans = 10                          # r_trans = amount transmitted
abx.s = 10                       # abxr_killr = amount of r killed by broad spectrum abx r
abx.r = 10                       # abxs_kills = amount of s killed by narrow spectrum abx s

