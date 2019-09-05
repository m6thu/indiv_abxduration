source('msm_util_rtnorm.R')
source('los_abx_matrix.R')

colo.table <- function(patient.matrix, los.array, total_prop, capacity_prop, prop_R, K){
    
    n.day = nrow(patient.matrix)
    n.bed = ncol(patient.matrix)
    
    number_of_patients = dim(los.array)[2]
    
    #total capacity for enterobacteriaceae growth (log)
    total_capacity= log(capacity_prop*exp(K)) 
    
    #existing population 
    total_existing= log(total_prop*exp(total_capacity)) #total number of Enterobacteriaceae is a proportion of the total capacity (log)
    
    total_bact = rnorm(number_of_patients, total_existing) #total number of Enterobacteriaceae for each patient (log)
    
    prop_R = rexp(1, 1/prop_R) #proportion of total Enterobacteriaceae that is resistant
    if (prop_R > 1) { prop_R = 1 }
    r_bact = log(prop_R*exp(total_bact)) #total number of resistant Enterobacteriaceae for each patient (log)
    s_bact = log(exp(total_bact) - exp(r_bact)) #total number of sensitive Enterobacteriaceae for each patient (log)
    
    S_Bactlevelstart = matrix(NA, n.day, n.bed)
    R_Bactlevelstart = matrix(NA, n.day, n.bed)
    
    # pad with NAs
    end_idx = 1
    for(i in 1:number_of_patients){
        S_Bactlevelstart[end_idx:(end_idx + los.array[2, i] - 1)] = c(s_bact[i], rep(NA, los.array[2, i]-1))
        R_Bactlevelstart[end_idx:(end_idx + los.array[2, i] - 1)] = c(r_bact[i], rep(NA, los.array[2, i]-1))
        end_idx = end_idx + los.array[2, i]
    }
    
    return(list(S_Bactlevelstart, R_Bactlevelstart)) # in log
}

# Update values for every day (define function)
nextDay <- function(patient.matrix, los.array, abx.matrix, colo.matrix, 
                    pi_ssr, total_prop, capacity_prop, K, r_thres, r_growth, r_trans, s_growth,
                    abx.s, abx.r, timestep){
    
    # K: loading capacity, and r_thres:threshold for infectiousness does not change with time
    pi_ssr = 1-(1-pi_ssr)^(1/timestep)
    
    #total capacity for enterobacteriaceae growth (log)
    total_capacity= log(capacity_prop*exp(K)) 
    
    #existing population 
    total_existing= log(total_prop*exp(total_capacity))
    
    S_table = colo.matrix[[1]] #in log
    R_table = colo.matrix[[2]] #in log
    
    # For each day (first day should be filled)
    for(i in 2:nrow(patient.matrix)){
        # calculate how many people has R above transmission level (log)
        r_num = sum(R_table[i-1,] > r_thres)
        # from number of r, calculate probability of transmission
        prop_r = 1-((1-pi_ssr)^r_num) 
        
        ###### Convert all log scale parameters into normal scale for addition, then convert back to log
        #for each person:
        for(j in 1:ncol(patient.matrix)){
            if(is.na(R_table[i, j])){ # pick any; S and R should be filled in same slots
                
                # roll for transmission
                roll = runif(1, 0, 1)
                # calculate effect of R logistic bacteria growth (abs)
                R_grow = r_growth*exp(R_table[i-1, j])*(1 - ((exp(R_table[i-1, j]) + exp(S_table[i-1, j]))/exp(total_capacity))) 
                # add effect of transmission if roll pass prob check and if previous R level is 0 (abs)
                R_trans = exp(r_trans)*(roll < prop_r)# & !R_table[i-1, j])
                # add effect of abx death if abx.matrix is r abx (== 2) (abs)
                R_abx = -(abx.matrix[i-1, j] == 2)*abx.r*R_grow
                # apply effects to current table (abs first because log of a negative number is NaN)
                R_table[i, j] = exp(R_table[i-1, j]) + R_grow + R_trans + R_abx
                # trim if value goes beyond range
                if(R_table[i, j] > exp(total_capacity)){
                    R_table[i, j] = exp(total_capacity) #abs
                }
                if(R_table[i, j] < 0){
                    R_table[i, j] = 0 #abs
                }
                R_table[i, j] = log(R_table[i, j]) #log
            }
            
            if(is.na(S_table[i, j])){ 
                # calculate effect of S logistic bacteria growth (in absolute numbers)
                S_grow = s_growth*exp(S_table[i-1, j])*(1 - ((exp(R_table[i-1, j]) + exp(S_table[i-1, j]))/exp(total_capacity)))
                # calculate effect of death antibiotics R and effect of death by abx S (abs)
                S_abx_s = -(abx.matrix[i-1, j] == 1)*abx.s*S_grow
                S_abx_r = -(abx.matrix[i-1, j] == 2)*abx.r*S_grow
                # apply effects (abs first because log of a negative number is NaN)
                S_table[i, j] = exp(S_table[i-1, j]) + S_grow + S_abx_s + S_abx_r
                # trim range
                if(S_table[i, j] > exp(total_capacity)){
                    S_table[i, j] = exp(total_capacity)
                }
                if(S_table[i, j] < 0){
                    S_table[i, j] = 0
                }
                S_table[i, j] = log(S_table[i, j]) #log
            }
        }
    }
    
    return(list(S_table, R_table))
}

diff_prevalence <- function(n.bed, max.los, p.infect, cum.r.1, p.r.day1,
                            K, total_prop, capacity_prop, prop_R, pi_ssr, r_thres, r_growth, r_trans, s_growth,
                            abx.s, abx.r, short_dur,long_dur){
    
    old = Sys.time() # get start time
    # DEBUG
    print(paste(n.bed, max.los, p.infect, cum.r.1, p.r.day1,
                K, total_prop, capacity_prop, prop_R, pi_ssr, r_thres, r_growth, r_trans, s_growth,
                abx.s, abx.r, short_dur,long_dur))
    
    timestep = 1
    n.day = 100
    iterations = 125
    
    #iter_totalR.no = matrix(NA, nrow = n.day, ncol = iterations)
    iter_totalR.thres = matrix(NA, nrow = n.day, ncol = iterations)
    
    for(iter in 1:iterations){
        
        matrixes = los.abx.table(n.bed=n.bed, n.day=n.day, max.los=max.los, 
                                 p.infect=p.infect, p.r.day1=p.r.day1, cum.r.1=cum.r.1, 
                                 meanDur= short_dur, timestep=timestep)
        patient.matrix=matrixes[[1]]
        abx.matrix=matrixes[[2]]
        los.array = summary.los(patient.matrix=patient.matrix)
        colo.matrix = colo.table(patient.matrix=patient.matrix, los.array=los.array, total_prop=total_prop, capacity_prop=capacity_prop, prop_R=prop_R, K=K)
        
        colo.matrix_filled_iter = nextDay(patient.matrix=patient.matrix, los.array=los.array, abx.matrix=abx.matrix, colo.matrix=colo.matrix, 
                                          pi_ssr=pi_ssr, total_prop=total_prop, capacity_prop=capacity_prop, K=K, r_thres=r_thres, r_growth=r_growth, r_trans=r_trans, s_growth=s_growth,
                                          abx.s=abx.s, abx.r=abx.r, timestep=timestep)
        # Summary - timestep by bed in absolute numbers
        df.R = data.frame(colo.matrix_filled_iter[[2]])
        
        # for number of people who reached R threshold on a day
        ##   sum of number of people per timestep that reach threshold 
        ##   make a matrix of sum of people per day (days by timestep)
        ##   daily means of number of people who reached R threshold 
        iter_totalR.thres[, iter]=rowMeans(matrix(rowSums(df.R >= r_thres), ncol=timestep, byrow=T))
        #print("end iteration loop")
    }
    totalR_thres_short = mean(rowSums(iter_totalR.thres[51:nrow(iter_totalR.thres),, drop=FALSE])/iterations/n.bed)
    
    #iter_totalR.no = matrix(NA, nrow = n.day, ncol = iterations)
    iter_totalR.thres = matrix(NA, nrow = n.day, ncol = iterations)
    
    for(iter in 1:iterations){
        
        matrixes = los.abx.table(n.bed=n.bed, n.day=n.day, max.los=max.los, 
                                 p.infect=p.infect, p.r.day1=p.r.day1, cum.r.1=cum.r.1, 
                                 meanDur= long_dur, timestep=timestep)
        patient.matrix=matrixes[[1]]
        abx.matrix=matrixes[[2]]
        los.array = summary.los(patient.matrix=patient.matrix)
        colo.matrix = colo.table(patient.matrix=patient.matrix, los.array=los.array, total_prop=total_prop, capacity_prop=capacity_prop, prop_R=prop_R,K=K)
        
        colo.matrix_filled_iter = nextDay(patient.matrix=patient.matrix, los.array=los.array, abx.matrix=abx.matrix, colo.matrix=colo.matrix, 
                                          pi_ssr=pi_ssr, total_prop=total_prop, capacity_prop=capacity_prop, K=K, r_thres=r_thres, r_growth=r_growth, r_trans=r_trans, s_growth=s_growth,
                                          abx.s=abx.s, abx.r=abx.r, timestep=timestep)
        
        #Summary 
        #for total units of R bacteria on a day
        df.R = data.frame(colo.matrix_filled_iter[[2]])
        #iter_totalR.no[, iter] = rowMeans(matrix(rowSums(df.R), ncol=timestep, byrow=T))
        #for number of people who reached R threshold on a day
        iter_totalR.thres[, iter] = rowMeans(matrix(rowSums(df.R >= r_thres), ncol=timestep, byrow=T))
        #print("end iteration loop")
    }
    #totalR_no_long = mean(rowSums(iter_totalR.no[51:nrow(iter_totalR.no),, drop=FALSE])/iterations/n.bed)
    totalR_thres_long = mean(rowSums(iter_totalR.thres[51:nrow(iter_totalR.thres),, drop=FALSE])/iterations/n.bed)
    
    # print elapsed time
    new = Sys.time() - old # calculate difference
    print(new) # print in nice format
    
    return(array(c(totalR_thres_long, totalR_thres_short, totalR_thres_long-totalR_thres_short)))
}

prevalence <- function(n.bed, max.los, p.infect, cum.r.1, p.r.day1,
                            K, total_prop, capacity_prop, prop_R, pi_ssr, r_thres, r_growth, r_trans, s_growth,
                            abx.s, abx.r, meanDur){
    
    old = Sys.time() # get start time
    # DEBUG
    print(paste(n.bed, max.los, p.infect, cum.r.1, p.r.day1,
                K, total_prop, capacity_prop, prop_R, pi_ssr, r_thres, r_growth, r_trans, s_growth,
                abx.s, abx.r, meanDur))
    
    timestep = 1
    n.day = 100
    iterations = 125
    
    iter_totalR.thres = matrix(NA, nrow = n.day, ncol = iterations)
    
    for(iter in 1:iterations){
        
        matrixes = los.abx.table(n.bed=n.bed, n.day=n.day, max.los=max.los, 
                                 p.infect=p.infect, p.r.day1=p.r.day1, cum.r.1=cum.r.1, 
                                 meanDur= meanDur, timestep=timestep)
        patient.matrix=matrixes[[1]]
        abx.matrix=matrixes[[2]]
        los.array = summary.los(patient.matrix=patient.matrix)
        colo.matrix = colo.table(patient.matrix=patient.matrix, los.array=los.array, total_prop=total_prop, capacity_prop=capacity_prop, prop_R=prop_R,K=K)
        
        colo.matrix_filled_iter = nextDay(patient.matrix=patient.matrix, los.array=los.array, abx.matrix=abx.matrix, colo.matrix=colo.matrix, 
                                          pi_ssr=pi_ssr, total_prop=total_prop, capacity_prop=capacity_prop, K=K, r_thres=r_thres, r_growth=r_growth, r_trans=r_trans, s_growth=s_growth,
                                          abx.s=abx.s, abx.r=abx.r, timestep=timestep)
        # Summary
        df.R = data.frame(colo.matrix_filled_iter[[2]])
        
        #for number of people who reached R threshold on a day
        iter_totalR.thres[, iter]=rowMeans(matrix(rowSums(df.R >= r_thres), ncol=timestep, byrow=T))
        #print("end iteration loop")
    }

    totalR_thres = mean(rowSums(iter_totalR.thres[51:nrow(iter_totalR.thres),, drop=FALSE])/iterations/n.bed)
    
    # print elapsed time
    new = Sys.time() - old # calculate difference
    print(new) # print in nice format
    
    return(totalR_thres)
}

parameters_prevalence_freq <- c("n.bed", "max.los", "p.infect", "cum.r.1", "p.r.day1", 
                     "K", "total_prop", "capacity_prop", "prop_R",
                     "pi_ssr", "r_thres", "r_growth", "r_trans", 's_growth',
                     "abx.s", "abx.r",
                     "meanDur")

parameters_diff_prevalence_freq <- c("n.bed", "max.los", "p.infect", "cum.r.1", "p.r.day1", 
                     "K", "total_prop", "capacity_prop", "prop_R",
                     "pi_ssr", "r_thres", "r_growth", "r_trans", 's_growth',
                     "abx.s", "abx.r",
                     "short_dur", "long_dur")

