setwd('/Users/moyin/Documents/nBox/git_projects/indiv_abxduration/')
rm(list=ls()) # Clean working environment

model='frequency'
source('default_params.R')
source(paste0('model_',model,'.R'))

#simple
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

n.bed =30
max.los = 10
p.infect=0.2
cum.r.1=791
p.r.day1=0.498
K=23.7
total_prop=0.33
prop_R=0.037
pi_ssr=0.059
r_trans=3.236
r_growth=0.0431
r_thres=8.1
s_growth=0.0059
abx.s=11.252
abx.r=0.00000001
short_dur=3.664
long_dur=16.9

matrixes=los.abx.table(n.bed, n.day, max.los, p.infect, p.r.day1, cum.r.1,meanDur, timestep=1)
patient.matrix=matrixes[[1]]
abx.matrix=matrixes[[2]]
los.array=summary.los(patient.matrix)
colo.matrix=colo.table(patient.matrix, los.array, total_prop, prop_R, r_thres, K)
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
