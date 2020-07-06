setwd('/Users/moyin/Documents/nBox/git_projects/indiv_abxduration/')
rm(list=ls()) # Clean working environment

model='frequency'
source(paste0('model_',model,'.R'))

#simple
all.filled = c()
for (t in 1:250) {
  parameters = c(runif(1, min=5, max=50),              #"n.bed", number of beds in the ward
                 runif(1, min=3, max=20),       #"mean.max.los", mean of length of stay
                 runif(1, min=0.01, max=0.8),    #"prob_StartBact_R",probability of initial carriage of resistant organisms
                 runif(1, min=0, max=1),         #"prop_S", proportion of S in the population of S and ss
                 runif(1, min=0, max=1),                  #"bif", bacterial interference factor
                 runif(1, min=0.001, max=0.3),            # "pi_ssr" probability of being transmitted r to ss (ss—> ssr)
                 runif(1, min=0.005, max=0.05),    # "repop.s1" probability of ss repopulated to S (Palleja, Nature Biology, 2018 on gut recovery ~9 months)
                 runif(1, min=0.002, max=0.02),         # "mu_r", probability of decolonisation (Haggai Bar-Yoseph, JAC, 2016, decreasing colonization rates from 76.7% (95% CI=69.3%–82.8%) at 1 month to 35.2% (95% CI=28.2%–42.9%) at 12 months of follow-up)
                 runif(1, min=0.1, max=0.5),          # "abx.s", probability of S becoming ss after being on narrow spectrum antibiotics
                 runif(1, min=0, max=0.00000001),            # "abx.r", probability of R becoming ss after being on broad spectrum antibiotics
                 runif(1, min=0.1, max=1),          # "p.infect", probability of being prescribed antibiotics
                 runif(1, min=30, max=300),         # admission day when cummulative prabability of HAI requiring abx.r is 1
                 runif(1, min=0.1, max=1),          # probability of being prescribed broad spectrum antibiotic on admission 
                 runif(1, min=3, max=7),           # "short_dur", mean short duration of antibiotics (normal distribution)
                 runif(1, min=14, max=21)          # "long_dur", mean long duration of antibiotics (normal distribution)
  )  
  parameters = c(12.83,11.194,0.73838,0.81,        
                 0.59,0.115218,0.02129,0.004052,0.1232,      
                 0.4112,0.6778,166.62,1,6.048,    
                 17.934)
  all.filled=diff_prevalence(n.bed = parameters[1], max.los = parameters[2], 
                                  prop_R=parameters[3], prop_S=parameters[4], 
                                  bif=parameters[5], pi_ssr=parameters[6], repop.s=parameters[7], mu=parameters[8], abx.s=parameters[9], abx.r=parameters[10],
                                  p.infect=parameters[11], cum.r.1=parameters[12], p.r.day1=parameters[13], short_dur=parameters[14], long_dur=parameters[15])
}

source('default_params.R')
matrixes=los.abx.table(n.bed, n.day, max.los,p.infect, p.r.day1, cum.r.1,meanDur, timestep=1)
patient.matrix=matrixes[[1]]
abx.matrix=matrixes[[2]]
los.array=summary.los(patient.matrix)
colo.matrix=colo.table(patient.matrix, los.array, prop_R, prop_S)

a=diff_prevalence(n.bed,max.los,prop_R, prop_S, bif, pi_ssr, repop.s, mu, abx.s, abx.r,
                  p.infect, cum.r.1, p.r.day1, short_dur, long_dur)[[3]]

prevalence(n.bed,max.los, prop_R, prop_S, bif, pi_ssr, repop.s, mu, abx.s, abx.r, p.infect, cum.r.1, p.r.day1, meanDur)

#binary 
matrixes=los.abx.table(n.bed, n.day, max.los,p.infect, p.r.day1, cum.r.1,meanDur, timestep=1)
patient.matrix=matrixes[[1]]
abx.matrix=matrixes[[2]]
los.array=summary.los(patient.matrix)
colo.matrix=colo.table(patient.matrix, los.array, prop_R, prop_r, prop_Sr, prop_S)
nextDay(patient.matrix, abx.matrix, colo.matrix, 
        pi_ssr, bif, mu, repop.r, repop.s, abx.r, abx.s, 
        timestep)

diff_prevalence(n.bed, max.los, prop_R, prop_r, prop_Sr, prop_S, bif, pi_ssr, repop.s, repop.r, mu, abx.s, abx.r, 
                p.infect, cum.r.1, p.r.day1, short_dur, long_dur)

prevalence(n.bed, max.los, prop_R, prop_r, prop_Sr, prop_S, bif, pi_ssr, repop.s,  repop.r, mu, abx.s, abx.r, 
           p.infect, cum.r.1, p.r.day1, meanDur)

#frequency
n.bed =19.1666666666667
max.los = 13.478835978836
p.infect=0.383333333333333
cum.r.1=147.857142857143
p.r.day1=0.65952380952381
K=15
total_prop=0.1
prop_R=0.8
pi_ssr=0.221690476190476
r_trans=7
r_growth=0.753703703703704
r_thres=  10
s_growth= 0.132142857142857
abx.s=  0.797354497354497
abx.r= 9.78835978835979e-15
short_dur=6.43915343915344
long_dur= 17.2037037037037

matrixes=los.abx.table(n.bed, n.day, max.los, p.infect, p.r.day1, cum.r.1,meanDur, timestep=1)
patient.matrix=matrixes[[1]]
abx.matrix=matrixes[[2]]
los.array=summary.los(patient.matrix)
colo.matrix=colo.table(patient.matrix, los.array, total_prop, prop_R, r_thres, K)
colo.matrix[[1]][1,]
S_table = colo.matrix[[1]] #in log
R_table = colo.matrix[[2]] #in log
n=nextDay(patient.matrix, los.array, abx.matrix, colo.matrix, 
          pi_ssr, total_prop, K, r_growth, r_thres, r_trans, s_growth,
          abx.s, abx.r, timestep=1)

diff_prevalence(n.bed, max.los, p.infect, cum.r.1, p.r.day1, K, total_prop, prop_R, pi_ssr, 
                r_trans, r_growth, r_thres, s_growth,
                abx.s, abx.r, short_dur,long_dur)

prevalence(n.bed, max.los, p.infect, cum.r.1, p.r.day1,
           K, total_prop, prop_R, pi_ssr, r_trans, r_growth,r_thres, s_growth,
           abx.s, abx.r, meanDur)
