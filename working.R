n.bed=35.6388888888889
mean.max.los=11.25 
p.infect=0.344444444444444
cum.r.1=5837.5
p.r.day1=0.744444444444444
K=19.5
total_prop=0.8
r_prop=0.8
pi_ssr=0.00416666666666667
r_thres=9.97222222222222
r_growth=1.05277777777778
s_growth=.05277777777778
r_trans=3.75
abx.s=5.75
abx.r=9.97222222222222
short_dur=6.88888888888889
long_dur=16.1388888888889

diff_prevalence(n.bed, mean.max.los, p.infect, cum.r.1, p.r.day1,
                            K, total_prop, r_prop, pi_ssr, r_thres, r_growth, r_trans, s_growth,
                            abx.s, abx.r, short_dur,long_dur)

matrixes=los.abx.table(n.bed, n.day=350, mean.max.los,
                          p.infect, p.r.day1, cum.r.1, 
                          meanDur=long_dur, timestep=1)
patient.matrix=matrixes[[1]]
abx.matrix=matrixes[[2]]
los.array=summary.los(patient.matrix)

colo.matrix=colo.table(patient.matrix, los.array, total_prop , r_prop,K)
nextDay(patient.matrix, los.array, abx.matrix, colo.matrix, 
                    pi_ssr, K, r_thres, r_growth, r_trans, s_growth,
                    abx.s, abx.r, timestep=1)
