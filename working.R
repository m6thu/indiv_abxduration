model='frequency'
source('default_params.R')
source(paste0('model_',model,'.R'))

#simple
diff_prevalence(n.bed,max.los, 
                            prop_R, prop_S_nonR, 
                            bif, pi_ssr, repop.s1, mu_r, abx.s, abx.r,
                            p.infect, cum.r.1, p.r.day1, short_dur, long_dur)

prevalence(n.bed,max.los, 
                prop_R, prop_S_nonR, 
                bif, pi_ssr, repop.s1, mu_r, abx.s, abx.r,
                p.infect, cum.r.1, p.r.day1, meanDur)

#binary 
diff_prevalence(n.bed, max.los, 
                            prop_R, prop_S_nonR, prop_Sr_inR, prop_sr_inR,
                            bif, pi_ssr, repop.s1, repop.s2, repop.r1, repop.r2,
                            mu1, mu2, mu_r, abx.s, abx.r, 
                            p.infect, cum.r.1, p.r.day1, short_dur, long_dur)

prevalence(n.bed, max.los, 
                prop_R, prop_S_nonR, prop_Sr_inR, prop_sr_inR,
                bif, pi_ssr, repop.s1, repop.s2, repop.r1, repop.r2,
                mu1, mu2, mu_r, abx.s, abx.r, 
                p.infect, cum.r.1, p.r.day1, meanDur)

#frequency
n.day=100
diff_prevalence(n.bed, max.los, p.infect, cum.r.1, p.r.day1,
                            K, total_prop, prop_R, pi_ssr, r_thres, r_growth, r_trans, s_growth,
                            abx.s, abx.r, short_dur,long_dur)

prevalence(n.bed, max.los, p.infect, cum.r.1, p.r.day1,
                K, total_prop, prop_R, pi_ssr, r_thres, r_growth, r_trans, s_growth,
                abx.s, abx.r, meanDur)
