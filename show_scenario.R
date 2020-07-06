##########################################################################
#######Effect of antibiotic duration on hospitalised patients############
#####################Show one run per scenario ##########################
#########################################################################
rm(list=ls()) # Clean working environment

#model can be "simple", "binary", or "frequency"
model <- "binary"
source(paste0("model_", model,".R"))

source('los_abx_matrix.R')
source('plot_functions/plot_scenario.R')

#adjust default variables here
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
                   cpe = list(n.bed = 10,       # n.bed= number of beds in the ward
                              max.los = 5,     # mean.max.los= mean of max length of stay (exponential distribution)
                              n.day = 30,
                              short_dur = 5,
                              long_dur = 15,
                              meanDur = 7,
                              prop_R = 0.2,     # Probability of being colonized with resistant strain on admission
                              pi_ssr = 0.05,     # pi_ssr= probability of R transmitting 
                              # simple cum.r.1 = 600, 
                              cum.r.1 = 300, 
                              p.infect = 0.3,   # p=probability of receiving antibiotic
                              p.r.day1 = 0.5, 
                              timestep = 1))

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
                              mu = 0.2,        # mu= probability of clearance of Sr to become S
                              repop.s = 0.1, 
                              repop.r= 0.3,    # probability of repopulation of Sr to become sR 
                              abx.s = 0.5,
                              abx.r = 0))

plot_scenario(model = model, scenario = 'cpe')

ggsave(paste0('../../../../Desktop/scenario_', model,'.jpeg'), units = 'cm', width = 32, height=25)
