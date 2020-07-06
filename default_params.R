############## Default parameters ##################

#common parameters 
common.para = list(esbl = list(n.bed = 20,       # n.bed= number of beds in the ward
                               max.los = 10,     # mean.max.los= mean of max length of stay (exponential distribution)
                               n.day = 300,
                               short_dur = 5,
                               long_dur = 15,
                               meanDur = 7,
                               prop_R = 0.8,     # Probability of being colonized with resistant strain on admission
                               pi_ssr = 0.2,     # pi_ssr= probability of R transmitting 
                               #simple cum.r.1 = 600, 
                               cum.r.1 = 900, 
                               p.infect = 0.5,   # p=probability of receiving antibiotic
                               p.r.day1 = 0.4, 
                               timestep = 1), 
                   cpe = list(n.bed = 20,       # n.bed= number of beds in the ward
                              max.los = 10,     # mean.max.los= mean of max length of stay (exponential distribution)
                              n.day = 300,
                              short_dur = 5,
                              long_dur = 15,
                              meanDur = 7,
                              prop_R = 0.8,     # Probability of being colonized with resistant strain on admission
                              pi_ssr = 0.2,     # pi_ssr= probability of R transmitting 
                              # simple cum.r.1 = 600, 
                              cum.r.1 = 900, 
                              p.infect = 0.5,   # p=probability of receiving antibiotic
                              p.r.day1 = 0.4, 
                              timestep = 1))

############# Simple model #################
simple.para = list(esbl = list(prop_S = 0.3,     # Proportion of large S within non-resistant states (S+s)
                               bif = 0.8,        # bacterial interference factor 
                               mu = 0.02,        # mu= probability of clearance of Sr to become S
                               repop.s = 0.4, 
                               abx.s = 0.3,
                               abx.r = 0.3), 
                   cpe = list(prop_S = 0.3,     # Proportion of large S within non-resistant states (S+s)
                              bif = 0.8,        # bacterial interference factor 
                              mu = 0.02,       # mu= probability of clearance of Sr to become S
                              repop.s = 0.2, 
                              abx.s = 0.3,
                              abx.r = 0))

############## Binary model ###################
binary.para = list(esbl = list(prop_S = 0.5,     # Proportion of large S within non-resistant states (S+s)
                               prop_Sr = 0.4,                  
                               prop_r = 0.5,    
                               bif = 0.7,        # bacterial interference factor
                               mu = 0.01,        # mu= probability of clearance of Sr to become S
                               repop.s = 0.01, 
                               repop.r= 0.03,    # probability of repopulation of Sr to become sR 
                               abx.s = 0.3,
                               abx.r = 0.3),
                   cpe = list(prop_S = 0.5,      # Proportion of large S within non-resistant states (S+s)
                              prop_Sr = 0.4,                  
                              prop_r = 0.5,    
                              bif = 0.7,        # bacterial interference factor
                              mu = 0.01,        # mu= probability of clearance of Sr to become S
                              repop.s = 0.01, 
                              repop.r= 0.03,    # probability of repopulation of Sr to become sR 
                              abx.s = 0.3,
                              abx.r = 0))

############# Frequency model #################
frequency.para = list(esbl = list(K = 22,           # gut holding capacity
                                  total_prop = 0.3, # mean of total starting amount of gut bacteria on log scale
                                  r_growth = 0.8,   # r_growth = growth constant for logistic growth
                                  s_growth = 0.1,
                                  r_trans = 4,      # r_thres = R threshold level for tranmissibility
                                  r_thres= 8, 
                                  abx.s = 0.7,
                                  abx.r = 0.7), 
                      cpe = list(K = 22,           # gut holding capacity
                                 total_prop = 0.3, # mean of total starting amount of gut bacteria on log scale
                                 r_growth = 0.8,   # r_growth = growth constant for logistic growth
                                 s_growth = 0.1,
                                 r_trans = 4,      # r_thres = R threshold level for tranmissibility
                                 r_thres= 8, 
                                 abx.s = 0.7,
                                 abx.r = 0))

if (model == 'simple') {
  
  para = list(common.para, simple.para)
  parameters_diff_prevalence = c("n.bed", "max.los", 
                                 "prop_R", "prop_S", 
                                 "bif", "pi_ssr", "repop.s", "mu", 
                                 "abx.s", "abx.r", "p.infect", "cum.r.1", 'p.r.day1', "short_dur", "long_dur")
  
} else if (model == 'binary') {
  
  para = list(common.para, binary.para)
  parameters_diff_prevalence = c("n.bed", "max.los", 
                                 "prop_R", "prop_r", "prop_Sr", "prop_S",
                                 "bif", "pi_ssr", "repop.s","repop.r",
                                 "mu", "abx.s", "abx.r",
                                 "p.infect", "cum.r.1", "p.r.day1",
                                 "short_dur", "long_dur")
  
} else if (model == 'frequency') {
  
  para = list(common.para, frequency.para)
  parameters_diff_prevalence = c("n.bed", "max.los", "p.infect", "cum.r.1", "p.r.day1", 
                                 "K", "total_prop", "prop_R",
                                 "pi_ssr", "r_trans", "r_growth", 'r_thres','s_growth',
                                 "abx.s", "abx.r",
                                 "short_dur", "long_dur")
}

#save each parameter in the list as an object
#for(i in 1:length(para)) assign(names(para)[i], para[[i]])
