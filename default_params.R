############## Default parameters ##################

#common parameters 
common.para = list(n.bed=30,        # n.bed= number of beds in the ward
                   max.los= 7,      # mean.max.los= mean of max length of stay (exponential distribution)
                   n.day= 365,
                   short_dur = 7,
                   long_dur = 14,
                   meanDur= 10,
                   prop_R=0.3,      # Probability of being colonized with resistant strain on admission
                   pi_ssr = 0.01,   # pi_ssr= probability of R transmitting 
                   cum.r.1=500, 
                   p.infect = 0.4,  # p=probability of receiving antibiotic
                   p.r.day1 =0.2)

############# Simple model #################
simple.para = list(prop_S = 0.8,    #Proportion of large S within non-resistant states (S+s)
                   bif = 0.5,       # bacterial interference factor 
                   mu = 0.01,       # mu= probability of clearance of Sr to become S
                   repop.s = 0.01, 
                   abx.s = 0.4,
                   abx.r = 0.2)

############## Binary model ###################
binary.para = list(prop_S = 0.8,     # Proportion of large S within non-resistant states (S+s)
                   prop_Sr = 0.7,                  
                   prop_r= 0.1,    
                   bif = 0.5,        # bacterial interference factor
                   mu = 0.01,        # mu= probability of clearance of Sr to become S
                   repop.s = 0.01, 
                   repop.r= 0.025,    # probability of repopulation of Sr to become sR 
                   abx.s = 0.3,
                   abx.r = 0.3)

############# Frequency model #################
frequency.para = list(K = 22,          # gut holding capacity
                      total_prop = 0.5,# mean of total starting amount of gut bacteria on log scale
                      r_growth = 0.025, # r_growth = growth constant for logistic growth
                      s_growth = 0.01,
                      r_trans = 5,     # r_thres = R threshold level for tranmissibility
                      r_thres=0.2, 
                      abx.s = 10,
                      abx.r = 10)

if (model == 'simple') {
  
  para = c(common.para, simple.para)
  parameters_diff_prevalence = c("n.bed", "max.los", 
                                 "prop_R", "prop_S", 
                                 "bif", "pi_ssr", "repop.s", "mu", 
                                 "abx.s", "abx.r", "p.infect", "cum.r.1", 'p.r.day1', "short_dur", "long_dur")
  
} else if (model == 'binary') {
  
  para = c(common.para, binary.para)
  parameters_diff_prevalence = c("n.bed", "max.los", 
                                  "prop_R", "prop_r", "prop_Sr", "prop_S",
                                  "bif", "pi_ssr", "repop.s","repop.r",
                                  "mu", "abx.s", "abx.r",
                                  "p.infect", "cum.r.1", "p.r.day1",
                                  "short_dur", "long_dur")
  
} else if (model == 'frequency') {
  
  para = c(common.para, frequency.para)
  parameters_diff_prevalence = c("n.bed", "max.los", "p.infect", "cum.r.1", "p.r.day1", 
                                  "K", "total_prop", "prop_R",
                                  "pi_ssr", "r_trans", "r_growth", 'r_thres','s_growth',
                                  "abx.s", "abx.r",
                                  "short_dur", "long_dur")
}

#save each parameter in the list as an object
for(i in 1:length(para)) assign(names(para)[i], para[[i]])