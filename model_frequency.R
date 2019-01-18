source('msm_util_rtnorm.R')

# generate a table of number of days we want to observe (rows) -
# against number of beds in the ward (columns), filled in with patient id numbers
patient.table <- function(n.bed, n.day, mean.max.los, timestep=1){
    
    #generate patient id numbers, the maximum number of patients possible is number of bed multiple by
    #number of days. This is to ensure there are enough total number of patients generated to fill table 
    n.patient <- n.bed*n.day 
    
    #vectorise the patient id to be used for filling in the patient.matrix
    patient.id <- 1:n.patient
    
    all_los <- ceiling(rexp(n.patient, 1/(mean.max.los*timestep)))
    all_los[all_los > 5*mean.max.los*timestep] <- 1
    sum_los <- cumsum(all_los)
    
    #make up a matrix of number of days we want to observe (rows) -
    #against number of beds in the ward (columns)
    patient.matrix <- matrix(NA, nrow=ceiling(n.day*timestep), ncol=n.bed)
    idx <- 1
    for(j in 1:n.bed){
        los_idx <- suppressWarnings(max(which(sum_los < n.day*timestep))) #Suppress warning that it creates -Inf
        # Handle case where first patient stays the whole observation duration
        if(los_idx == -Inf){
            los_idx <- 1
            #print(idx:(idx+length(los)))
            patient.matrix[, j] <- rep(idx, n.day*timestep)
            #print('pat')
            #print(patient.matrix[, j])
            idx <- idx+1
            all_los <- all_los[-(1)]
        }else{
            los <- all_los[1:los_idx]
            #print(idx:(idx+length(los)))
            patient.matrix[, j] <- rep(idx:(idx+length(los)), c(los, n.day*timestep-sum(los)))
            #print('pat')
            #print(patient.matrix[, j])
            idx <- idx+length(los)+1
            all_los <- all_los[-(1:(los_idx+1))]
        }
        sum_los <- cumsum(all_los)
    }
    
    return(patient.matrix)
}

#frequency summary of patient.matrix - patient id against number of days of stay for each patient
summary.los <- function(patient.matrix){    
    
    # Summarize how often each patient.id is encountered to get days for each id
    los.dur <- table(patient.matrix)
    los_duration <- array(dim = c(2, length(los.dur)))
    # Attach patient ID on 1st row
    los_duration[1,] <- 1:length(los.dur)
    # Put summary of days on 2nd row
    los_duration[2,] <- los.dur
    
    return(los_duration)
}

abx.table <- function(patient.matrix, los.array, p.s, p.r.day1, p.r.dayafter,
                      meanDur.s, meanDur.r, sdDur, timestep=1){
    
    p.s <- p.s/timestep
    p.r.day1 <- p.r.day1/timestep
    p.r.dayafter <- p.r.dayafter/timestep
    # Check assumption that possibilities are, no abx, has s abx, or has r abx on first day
    stopifnot(p.s+p.r.day1 <= 1)
    
    #generate antibiotic use table
    #number of days of s antibiotic is randomly drawn from a truncated normal distribution
    abx_days.s <- round(rtnorm(ncol(los.array), mean=meanDur.s*timestep, sd=sdDur*timestep, lower=1))
    #number of days of r antibiotic is randomly drawn from a truncated normal distribution
    abx_days.r <- round(rtnorm(ncol(los.array), mean=meanDur.r*timestep, sd=sdDur*timestep, lower=1))
    #number of days of r antibiotic for days after drawn from tunced norm dist
    abx_r.after <- round(rtnorm(length(patient.matrix), mean=meanDur.r*timestep, sd=sdDur*timestep, lower=1))
    r_idx <- 1 # R indices start at 1
    # Unit test - check distribution of abx distribution
    # hist(abx_days.s, breaks=20)
    # Unit test - compare cases that will enter padding if-else
    # abx_days > los.array[2, ]
    
    # abx.matrix should be same size as patient.matrix
    abx.matrix <- matrix(NA, nrow=nrow(patient.matrix), ncol=ncol(patient.matrix))
    idx_end <- 1
    # For each patient
    for(i in 1:ncol(los.array)){
        
        # maxiumum number of days for a particular patient from los vector
        max_days <- los.array[2, i]
        # number of abx days for that particular patient from generated number
        abx_s <- abx_days.s[i]
        abx_r <- abx_days.r[i]
        
        # Initial treatment value derived from probability for p.r and p.s
        rand <- runif(1, 0, 1)
        if (rand < p.s){ # check if p.s starting
            if(abx_s > max_days){
                # if the number of generated abx is longer than the los of that person
                # have them take abx everyday for their stay
                abx.matrix[idx_end:(idx_end+max_days-1)] <- rep(1, max_days)
            }else{
                # else take abx for abx days and pad to fit max_days
                abx.matrix[idx_end:(idx_end+max_days-1)] <- c(rep(1, abx_s), rep(0, max_days-abx_s))
            }
        }else if(rand < (p.s+p.r.day1)){ # case for p.r starting, rand >= p.s already checked in previous if case
            if(abx_r > max_days){
                # if the number of generated abx is longer than the los of that person
                # have them take abx everyday for their stay
                abx.matrix[idx_end:(idx_end+max_days-1)] <- rep(2, max_days)
            }else{
                # else take abx for abx days and pad to fit max_days
                abx.matrix[idx_end:(idx_end+max_days-1)] <- c(rep(2, abx_r), rep(0, max_days-abx_r))
            }
        }else{
            # no abx taken for that person
            abx.matrix[idx_end:(idx_end+max_days-1)] <- rep(0, max_days)
        }
        
        #Every day has a chance of starting abx.r by a fixed probability
        roll.r <- runif(max_days, 0, 1)
        start.r <- roll.r < p.r.dayafter
        #print(sum(start.r))
        where.r <- which(start.r == 1)
        #print(paste('where.r', where.r))
        # Create an array of p.r starting based on this
        # Merge starting vector and random p.r start vector (giving priority to p.r random start)
        # for each location that rolled chance that abx.r starts
        if(sum(start.r)){
            for(j in 1:sum(start.r)){
                #print(j)
                # replace that location for the length of abx_r.after drawn from norm distribution
                start_idx <- idx_end+where.r[j]-1
                #print(paste(start_idx, abx_r.after[r_idx], length(abx_r.after), r_idx))
                end_idx <- start_idx+abx_r.after[r_idx]-1
                #print(paste(end_idx, idx_end, max_days))
                if(end_idx > (idx_end+max_days-1)){ # abx duration exceeds max_days for that person, replace the rest, escape loop
                    abx.matrix[start_idx:(idx_end+max_days-1)] <- rep(2, max_days-where.r[j]+1)
                    break
                }else{ # abx duration does not exceed max_days for that person
                    abx.matrix[start_idx:end_idx] <- rep(2, abx_r.after[r_idx])
                }
                # move to next abx_r.after drawn
                r_idx <- r_idx+1
            }
        }
        
        # move starting position to end of previous patient
        idx_end <- idx_end+max_days
    }
    
    return(abx.matrix)
}

# Defaults from Rene's data ini_16S_log
colo.table <- function(patient.matrix, los.array, t_mean, t_sd, r_mean, r_sd){
    
    n.day <- nrow(patient.matrix)
    n.bed <- ncol(patient.matrix)
    
    number_of_patients <- dim(los.array)[2]
    # Unit test example: mean(t_mean), sd(t_sd)
    total_bact <- rnorm(number_of_patients, t_mean, t_sd)
    r_bact <- rnorm(number_of_patients, r_mean, r_sd)
    
    # since both are on log scale, to sum them together would need logsumexp()
    s_bact <- exp(total_bact) - exp(r_bact)
    # spin until no minus values for log... probably not best way
    while(min(s_bact) < 0){
        spin_n <- sum(s_bact < 0)
        spin_total <- rnorm(spin_n, t_mean, t_sd)
        spin_r <- rnorm(spin_n, r_mean, r_sd)
        s_bact[s_bact < 0] <- exp(spin_total) - exp(spin_r)
    }
    # convert back to log, log part of logsumexp()
    s_bact <- log(s_bact)
    
    S_Bactlevelstart <- matrix(NA, n.day, n.bed)
    R_Bactlevelstart <- matrix(NA, n.day, n.bed)
    
    # pad with NAs
    end_idx <- 1
    for(i in 1:number_of_patients){
        S_Bactlevelstart[end_idx:(end_idx + los.array[2, i] - 1)] <- c(s_bact[i], rep(NA, los.array[2, i]-1))
        R_Bactlevelstart[end_idx:(end_idx + los.array[2, i] - 1)] <- c(r_bact[i], rep(NA, los.array[2, i]-1))
        end_idx = end_idx + los.array[2, i]
    }
    
    return(list(S_Bactlevelstart, R_Bactlevelstart))
}

# 4. Update values for every day (define function)
nextDay <- function(patient.matrix, los.array, abx.matrix, colo.matrix, 
                    pi_r, K, r_thres, r_growth, r_trans, 
                    abxr_killr, abxr_kills, abxs_kills, timestep=1){
    
    # K: loading capacity, and r_thres:threshold for infectiousness does not change with time
    pi_r <- pi_r/timestep
    r_growth <- r_growth/timestep
    r_trans <- r_trans/timestep
    abxr_killr <- abxr_killr/timestep
    abxr_kills <- abxr_kills/timestep
    
    S_table <- colo.matrix[[1]]
    R_table <- colo.matrix[[2]]
    
    # For each day (first day should be filled)
    for(i in 2:nrow(patient.matrix)){
        # calculate how many people has R above transmission level
        r_num <- sum(R_table[i-1,] > r_thres)
        # from number of r, calculate probability of transmission
        prob_r <- 1-((1-pi_r)^r_num)
        
        ###### Convert all log scale parameters into normal scale for addition, then convert back to log
        #for each person:
        for(j in 1:ncol(patient.matrix)){
            if(is.na(R_table[i, j])){ # pick any; S and R should be filled in same slots
                # roll for transmission
                roll <- runif(1, 0, 1)
                # calculate effect of R logistic bacteria growth 
                R_grow <- r_growth*exp(R_table[i-1, j])*(1 - (exp(R_table[i-1, j]) + exp(S_table[i-1, j]))/K)
                # add effect of transmission if roll pass prob check and if previous R level is 0
                R_trans <- exp(r_trans)*((roll < prob_r))# & !R_table[i-1, j])
                # add effect of abx death if abx.matrix is r abx (== 2)
                R_abx <- -(abx.matrix[i-1, j] == 2)*exp(abxr_killr)
                # apply effects to current table
                R_table[i, j] <- exp(R_table[i-1, j]) + R_grow + R_trans + R_abx
                # trim if value goes beyond range
                if(R_table[i, j] > K){
                    R_table[i, j] <- K
                }
                if(R_table[i, j] < 0){
                    R_table[i, j] <- 0
                }
                R_table[i, j] <- log(R_table[i, j])
                if(R_table[i, j] < 0){
                    R_table[i, j] <- 0
                }
                
                # calculate effect of S logistic bacteria growth 
                S_grow <- r_growth*exp(S_table[i-1, j])*(1 - (exp(R_table[i-1, j]) + exp(S_table[i-1, j]))/K)
                # calculate effect of death antibiotics R and effect of death by abx S
                S_abx_s <- -(abx.matrix[i-1, j] == 1)*exp(abxs_kills)
                S_abx_r <- -(abx.matrix[i-1, j] == 2)*exp(abxr_kills)
                # apply effects
                S_table[i, j] <- exp(S_table[i-1, j]) + R_grow + S_abx_s + S_abx_r
                # trim range
                if(S_table[i, j] > K){
                    S_table[i, j] <- K
                }
                if(S_table[i, j] < 0){
                    S_table[i, j] <- 0
                }
                S_table[i, j] <- log(S_table[i, j])
                if(S_table[i, j] < 0){
                    S_table[i, j] <- 0
                }
            }
        }
        
    }
    
    return(list(S_table, R_table))
}

diff_prevalence <- function(n.bed, mean.max.los, p.s, p.r.day1, p.r.dayafter,
                            K, t_mean, t_sd, r_mean, r_sd,
                            pi_r, r_thres, r_growth, r_trans, 
                            abxr_killr, abxr_kills, abxs_kills,
                            short_dur.s, long_dur.s, short_dur.r, long_dur.r, sdDur){
    
    old <- Sys.time() # get start time
    # DEBUG
    print(paste(n.bed, mean.max.los, p.s, p.r.day1, p.r.dayafter,
                K, t_mean, t_sd, r_mean, r_sd,
                pi_r, r_thres, r_growth, r_trans, 
                abxr_killr, abxr_kills, abxs_kills,
                short_dur.s, long_dur.s, short_dur.r, long_dur.r, sdDur))
    
    timestep <- 2
    n.day <- 500
    iterations <- 100
    
    iter_totalR.no <- matrix(NA, nrow = n.day*timestep, ncol = iterations)
    iter_totalR.thres <- matrix(NA, nrow = n.day*timestep, ncol = iterations)
    
    for(iter in 1:iterations){
        
        patient.matrix <- patient.table(n.bed, n.day, mean.max.los, timestep)
        los.array <- summary.los(patient.matrix)
        abx.matrix <- abx.table(patient.matrix, los.array, p.s=p.s, p.r.day1=p.r.day1, p.r.dayafter=p.r.dayafter,
                                meanDur.s=short_dur.s, meanDur.r=short_dur.r, sdDur=sdDur, timestep=timestep)
        colo.matrix <- colo.table(patient.matrix, los.array, t_mean, t_sd, r_mean, r_sd)
        
        colo.matrix_filled_iter <- nextDay(patient.matrix, los.array, abx.matrix, colo.matrix, 
                                           pi_r, K, r_thres, r_growth, r_trans, 
                                           abxr_killr, abxr_kills, abxs_kills, timestep)
        # Summary
        df.R <- data.frame(colo.matrix_filled_iter[[2]])
        print(df.R)
        iter_totalR.no[, iter] <- rowSums(df.R)
        
        #for number of people who reached R threshold on a day
        iter_totalR.thres[, iter]<-rowSums(df.R >= r_thres)
        #print("end iteration loop")
    }
    totalR_no_short <- mean(rowSums(iter_totalR.no)/iterations/n.bed)
    totalR_thres_short <- mean(rowSums(iter_totalR.thres)/iterations/n.bed)
    
    iter_totalR.no <- matrix(NA, nrow = n.day, ncol = iterations)
    iter_totalR.thres <- matrix(NA, nrow = n.day, ncol = iterations)
    
    for(iter in 1:iterations){
        
        patient.matrix <- patient.table(n.bed, n.day, mean.max.los, timestep)
        los.array <- summary.los(patient.matrix)
        abx.matrix <- abx.table(patient.matrix, los.array, p.s=p.s, p.r.day1=p.r.day1, p.r.dayafter=p.r.dayafter,
                                meanDur.s=short_dur.s, meanDur.r=short_dur.r, sdDur=sdDur, timestep=timestep)
        colo.matrix <- colo.table(patient.matrix, los.array, t_mean, t_sd, r_mean, r_sd)
        
        colo.matrix_filled_iter <- nextDay(patient.matrix, los.array, abx.matrix, colo.matrix, 
                                           pi_r, K, r_thres, r_growth, r_trans, 
                                           abxr_killr, abxr_kills, abxs_kills, timestep)

        #Summary 
        #for total units of R bacteria on a day
        df.R <- data.frame(colo.matrix_filled_iter[[2]])
        iter_totalR.no[, iter] <- rowSums(df.R)
        
        #for number of people who reached R threshold on a day
        iter_totalR.thres[, iter] <- rowSums(df.R >= r_thres)
        #print("end iteration loop")
    }
    totalR_no_long <- mean(rowSums(iter_totalR.no)/iterations/n.bed)
    totalR_thres_long <- mean(rowSums(iter_totalR.thres)/iterations/n.bed)
    
    # print elapsed time
    new <- Sys.time() - old # calculate difference
    print(new) # print in nice format
    
    return(array(c((totalR_no_long - totalR_no_short), (totalR_thres_long-totalR_thres_short))))
}
res.names <- c(paste("No R per bed"),paste("R Thres per bed"))

parameters_frequency <- c("n.bed", "mean.max.los", "p.s", "p.r.day1", "p.r.dayafter",
                          "K", "t_mean", "t_sd", "r_mean", "r_sd",
                          "pi_r", "r_thres", "r_growth", "r_trans",
                          "abxr_killr", "abxr_kills", "abxs_kills",
                          "short_dur.s", "long_dur.s", "short_dur.r", "long_dur.r", "sdDur")

