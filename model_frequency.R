source('msm_util_rtnorm.R')
source('los_abx_matrix.R')

# Defaults from Rene's data ini_16S_log
colo.table <- function(patient.matrix, los.array, t_mean, r_mean){

    n.day = nrow(patient.matrix)
    n.bed = ncol(patient.matrix)
    
    number_of_patients = dim(los.array)[2]
    # Unit test example: mean(t_mean), sd()
    total_bact = rnorm(number_of_patients, t_mean )
    r_bact = rnorm(number_of_patients, r_mean )
    
    # since both are on log scale, to sum them together would need logsumexp()
    s_bact = exp(total_bact) - exp(r_bact)
    # spin until no minus values for log... probably not best way
    while(min(s_bact) < 0){
        spin_n = sum(s_bact < 0)
        spin_total = rnorm(spin_n, t_mean)
        spin_r = rnorm(spin_n, r_mean)
        s_bact[s_bact < 0] = exp(spin_total) - exp(spin_r)
    }
    # convert back to log, log part of logsumexp()
    s_bact = log(s_bact)
    
    S_Bactlevelstart = matrix(NA, n.day, n.bed)
    R_Bactlevelstart = matrix(NA, n.day, n.bed)
    
    # pad with NAs
    end_idx = 1
    for(i in 1:number_of_patients){
        S_Bactlevelstart[end_idx:(end_idx + los.array[2, i] - 1)] = c(s_bact[i], rep(NA, los.array[2, i]-1))
        R_Bactlevelstart[end_idx:(end_idx + los.array[2, i] - 1)] = c(r_bact[i], rep(NA, los.array[2, i]-1))
        end_idx = end_idx + los.array[2, i]
    }
    
    return(list(S_Bactlevelstart, R_Bactlevelstart))
}

# Update values for every day (define function)
nextDay <- function(patient.matrix, los.array, abx.matrix, colo.matrix, 
                    pi_r, K, r_thres, r_growth, r_trans, 
                    abx.s, abx.r, timestep){
    
    # K: loading capacity, and r_thres:threshold for infectiousness does not change with time
    pi_r = 1-(1-pi_r)^(1/timestep)
    
    if (abx.r < 0.1) { 
        abxr_killr = abx.r
        abxr_kills = abx.s
        abxs_kills = abx.s
    } else {
        abxr_killr = abx.s
        abxr_kills = abx.s
        abxs_kills = abx.s
    }
    
    S_table = colo.matrix[[1]]
    R_table = colo.matrix[[2]]
    
    # For each day (first day should be filled)
    for(i in 2:nrow(patient.matrix)){
        # calculate how many people has R above transmission level
        r_num = sum(R_table[i-1,] > r_thres)
        # from number of r, calculate probability of transmission
        prob_r = 1-((1-pi_r)^r_num)
        
        ###### Convert all log scale parameters into normal scale for addition, then convert back to log
        #for each person:
        for(j in 1:ncol(patient.matrix)){
            if(is.na(R_table[i, j])){ # pick any; S and R should be filled in same slots
                # roll for transmission
                roll = runif(1, 0, 1)
                # calculate effect of R logistic bacteria growth 
                R_grow = r_growth*exp(R_table[i-1, j])*(1 - (exp(R_table[i-1, j]) + exp(S_table[i-1, j]))/K)
                # add effect of transmission if roll pass prob check and if previous R level is 0
                R_trans = exp(r_trans)*(roll < prob_r)# & !R_table[i-1, j])
                # add effect of abx death if abx.matrix is r abx (== 2)
                R_abx = -(abx.matrix[i-1, j] == 2)*exp(abxr_killr)
                # apply effects to current table
                R_table[i, j] = exp(R_table[i-1, j]) + R_grow + R_trans + R_abx
                # trim if value goes beyond range
                if(R_table[i, j] > K){
                    R_table[i, j] = K
                }
                if(R_table[i, j] < 0){
                    R_table[i, j] = 0
                }
                R_table[i, j] = log(R_table[i, j])
                if(R_table[i, j] < 0){
                    R_table[i, j] = 0
                }
                
                # calculate effect of S logistic bacteria growth 
                S_grow = r_growth*exp(S_table[i-1, j])*(1 - (exp(R_table[i-1, j]) + exp(S_table[i-1, j]))/K)
                # calculate effect of death antibiotics R and effect of death by abx S
                S_abx_s = -(abx.matrix[i-1, j] == 1)*exp(abxs_kills)
                S_abx_r = -(abx.matrix[i-1, j] == 2)*exp(abxr_kills)
                # apply effects
                S_table[i, j] = exp(S_table[i-1, j]) + R_grow + S_abx_s + S_abx_r
                # trim range
                if(S_table[i, j] > K){
                    S_table[i, j] = K
                }
                if(S_table[i, j] < 0){
                    S_table[i, j] = 0
                }
                S_table[i, j] = log(S_table[i, j])
                if(S_table[i, j] < 0){
                    S_table[i, j] = 0
                }
            }
        }
        
    }
    
    return(list(S_table, R_table))
}

diff_prevalence <- function(n.bed, mean.max.los, p.infect, cum.r.1, p.r.day1,
                            K, t_mean, r_mean, pi_r, r_thres, r_growth, r_trans, 
                            abx.s, abx.r, short_dur,long_dur){
    
    old = Sys.time() # get start time
    # DEBUG
    print(paste(n.bed, mean.max.los, p.infect, cum.r.1, p.r.day1,
                K, t_mean, r_mean, pi_r, r_thres, r_growth, r_trans, 
                abx.s, abx.r, short_dur,long_dur))
    
    timestep = 1
    n.day = 350
    iterations = 10
    
    iter_totalR.no = matrix(NA, nrow = n.day*timestep, ncol = iterations)
    iter_totalR.thres = matrix(NA, nrow = n.day*timestep, ncol = iterations)
    
    for(iter in 1:iterations){
        
        matrixes = los.abx.table(n.bed=n.bed, n.day=n.day, mean.max.los=mean.max.los, 
                                 p.infect=p.infect, p.r.day1=p.r.day1, cum.r.1=cum.r.1, 
                                 meanDur= short_dur, timestep=timestep)
        patient.matrix=matrixes[[1]]
        abx.matrix=matrixes[[2]]
        los.array = summary.los(patient.matrix=patient.matrix)
        colo.matrix = colo.table(patient.matrix=patient.matrix, los.array=los.array, t_mean=t_mean, r_mean=r_mean)
        
        colo.matrix_filled_iter = nextDay(patient.matrix=patient.matrix, los.array=los.array, abx.matrix=abx.matrix, colo.matrix=colo.matrix, 
                                           pi_r=pi_r, K=K, r_thres=r_thres, r_growth=r_growth, r_trans=r_trans, 
                                           abx.s=abx.s, abx.r=abx.r, timestep=timestep)
        # Summary
        df.R = data.frame(colo.matrix_filled_iter[[2]])
        #print(df.R)
        iter_totalR.no[, iter] = rowSums(df.R)
        
        #for number of people who reached R threshold on a day
        iter_totalR.thres[, iter]=rowSums(df.R >= r_thres)
        #print("end iteration loop")
    }
    totalR_no_short = mean(rowSums(iter_totalR.no)/iterations/n.bed)
    totalR_thres_short = mean(rowSums(iter_totalR.thres)/iterations/n.bed)
    
    iter_totalR.no = matrix(NA, nrow = n.day*timestep, ncol = iterations)
    iter_totalR.thres = matrix(NA, nrow = n.day*timestep, ncol = iterations)
    
    for(iter in 1:iterations){
        
        matrixes = los.abx.table(n.bed=n.bed, n.day=n.day, mean.max.los=mean.max.los, 
                                 p.infect=p.infect, p.r.day1=p.r.day1, cum.r.1=cum.r.1, 
                                 meanDur= long_dur, timestep=timestep)
        patient.matrix=matrixes[[1]]
        abx.matrix=matrixes[[2]]
        los.array = summary.los(patient.matrix=patient.matrix)
        colo.matrix = colo.table(patient.matrix=patient.matrix, los.array=los.array, t_mean=t_mean, r_mean=r_mean)
        
        colo.matrix_filled_iter = nextDay(patient.matrix=patient.matrix, los.array=los.array, abx.matrix=abx.matrix, colo.matrix=colo.matrix, 
                                           pi_r=pi_r, K=K, r_thres=r_thres, r_growth=r_growth, r_trans=r_trans, 
                                           abx.s=abx.s, abx.r=abx.r, timestep=timestep)

        #Summary 
        #for total units of R bacteria on a day
        df.R = data.frame(colo.matrix_filled_iter[[2]])
        iter_totalR.no[, iter] = rowSums(df.R)
        
        #for number of people who reached R threshold on a day
        iter_totalR.thres[, iter] = rowSums(df.R >= r_thres)
        #print("end iteration loop")
    }
    totalR_no_long = mean(rowSums(iter_totalR.no)/iterations/n.bed)
    totalR_thres_long = mean(rowSums(iter_totalR.thres)/iterations/n.bed)
    
    # print elapsed time
    new = Sys.time() - old # calculate difference
    print(new) # print in nice format
    
    return(array(c((totalR_no_long - totalR_no_short), (totalR_thres_long-totalR_thres_short))))
}
res.names <- c(paste("No R per bed"),paste("R Thres per bed"))

parameters_freq <- c("n.bed", "mean.max.los", "p.infect", "cum.r.1", "p.r.day1", 
                          "K", "t_mean", "r_mean",
                          "pi_r", "r_thres", "r_growth", "r_trans",
                          "abx.s", "abx.r",
                          "short_dur", "long_dur")

