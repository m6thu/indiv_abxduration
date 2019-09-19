############## Shared parameters ##################
#fixed parameters 
n.bed=30                        # n.bed= number of beds in the ward
max.los=7                       # mean.max.los= mean of max length of stay (exponential distribution)
n.day=365

############# Simple model #################
short_dur = 5
long_dur = 14
meanDur=10

p.infect = 0.4                   # p=probability of receiving antibiotic
p.r.day1 =0.3
prop_R = 0.3                     # Probability of being colonized with resistant strain on admission
prop_S_nonR = 0.8                # Proportion of large S within non-resistant states (S+s)
bif = 0.5                        # bacterial interference factor 
mu = 0.01                        # mu1= probability of clearance of Sr to become S
pi_ssr = 0.0002                  # pi_s= probability of R transmitting to ss
repop.s = 0.01
cum.r.1=500

############## Binary model ###################
prop_S_nonR = 0.8                # Proportion of large S within non-resistant states (S+s)
prop_Sr_inR = 0.7                  
prop_sr_inR = 0.1                     
repop.r= 0.01                    # probability of repopulation of Sr to become sR 

############# Frequency model #################
K = 25                           # gut holding capacity
total_prop = 0.4                 # mean of total starting amount of gut bacteria on log scale
#t_sd = 1.1218                   # sd of total starting amount of gut bacteria on log scale
#r_sd = 1.8921                   # sd of starting amount of resistant gut bacteria on log scale
r_growth = 0.02                  # r_growth = growth constant for logistic growth
s_growth = 0.02
r_mean = 0.12                    # r_thres = R threshold level for tranmissibility
r_thres=0.5

if (model=='simple' | model=='binary'){
    abx.s = 0.2                  # probability of clearing S to become ss under antibiotic treatment 
    abx.r = 0.2
} else {
    abx.s = 12                   # abxr_killr = amount of r killed by broad spectrum abx r
    abx.r = 12                   # abxs_kills = amount of s killed by narrow spectrum abx s
}
    
