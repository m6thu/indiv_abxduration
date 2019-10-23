############## Shared parameters ##################
#fixed parameters 
n.bed=30                        # n.bed= number of beds in the ward
max.los=7                       # mean.max.los= mean of max length of stay (exponential distribution)
n.day=365

############# Simple model #################
short_dur = 7
long_dur = 14
meanDur=10

p.infect = 0.4                   # p=probability of receiving antibiotic
p.r.day1 =0.2
prop_R = 0.3                     # Probability of being colonized with resistant strain on admission
prop_S = 0.8                # Proportion of large S within non-resistant states (S+s)
bif = 0.5                        # bacterial interference factor 
mu = 0.01                        # mu1= probability of clearance of Sr to become S
pi_ssr = 0.01                  # pi_s= probability of R transmitting to ss
repop.s = 0.01
cum.r.1=500



############## Binary model ###################
prop_S = 0.8                # Proportion of large S within non-resistant states (S+s)
prop_Sr = 0.7                  
prop_r= 0.1                     
repop.r= 0.01                    # probability of repopulation of Sr to become sR 

############# Frequency model #################
K = 20                           # gut holding capacity
total_prop = 0.4                 # mean of total starting amount of gut bacteria on log scale
#t_sd = 1.1218                   # sd of total starting amount of gut bacteria on log scale
#r_sd = 1.8921                   # sd of starting amount of resistant gut bacteria on log scale
r_growth = 0.05                  # r_growth = growth constant for logistic growth
s_growth = 0.015
r_trans = 5                    # r_thres = R threshold level for tranmissibility
r_thres=0.2
