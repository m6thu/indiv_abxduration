model='frequency'
source('default_params.R')
source(paste0('model_',model,'.R'))

#simple
diff_prevalence(n.bed,max.los,prop_R, prop_S_nonR, bif, pi_ssr, repop.s, mu, abx.s, abx.r,
                p.infect, cum.r.1, p.r.day1, short_dur, long_dur)

prevalence(n.bed,max.los, 
           prop_R, prop_S_nonR, 
           bif, pi_ssr, repop.s, mu, abx.s, abx.r,
           p.infect, cum.r.1, p.r.day1, meanDur)

#binary 
n.day=300
max.los=7
timestep=1
matrixes=los.abx.table(n.bed, n.day, max.los,
                       p.infect, p.r.day1, cum.r.1, 
                       meanDur, timestep)
patient.matrix=matrixes[[1]]
abx.matrix=matrixes[[2]]
los.array=summary.los(patient.matrix)
colo.matrix=colo.table(patient.matrix, los.array, prop_R, prop_S_nonR, prop_Sr_inR, prop_sr_inR)

diff_prevalence(n.bed, max.los, 
                prop_R, prop_S_nonR, prop_Sr_inR, prop_sr_inR,
                bif, pi_ssr, repop.s, repop.r, mu, abx.s, abx.r, 
                p.infect, cum.r.1, p.r.day1, short_dur, long_dur)

prevalence(n.bed, max.los, 
           prop_R, prop_S_nonR, prop_Sr_inR, prop_sr_inR,
           bif, pi_ssr, repop.s,  repop.r, mu, abx.s, abx.r, 
           p.infect, cum.r.1, p.r.day1, meanDur)

#frequency
n.day=300
max.los=7
timestep=1
abx.s=11
abx.r=11
matrixes=los.abx.table(n.bed, n.day, max.los,
                       p.infect, p.r.day1, cum.r.1, 
                       meanDur, timestep)
patient.matrix=matrixes[[1]]
abx.matrix=matrixes[[2]]
los.array=summary.los(patient.matrix)
colo.matrix=colo.table(patient.matrix, los.array, total_prop, prop_R, r_mean, K)
n=nextDay(patient.matrix, los.array, abx.matrix, colo.matrix, 
          pi_ssr, total_prop, K, r_mean, r_growth, r_thres, s_growth,
          abx.s, abx.r, timestep)

diff_prevalence(n.bed, max.los, p.infect, cum.r.1, p.r.day1,
                K, total_prop, prop_R, pi_ssr, 
                r_mean, r_growth, r_thres, s_growth,
                abx.s, abx.r, short_dur,long_dur)

prevalence(n.bed, max.los, p.infect, cum.r.1, p.r.day1,
           K, total_prop, prop_R, pi_ssr, r_mean, r_growth,r_thres, s_growth,
           abx.s, abx.r, meanDur)
