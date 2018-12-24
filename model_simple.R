require(msm)

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

abx.table <- function(patient.matrix, los.array, p, meanDur, sdDur, timestep=1){
    
    #generate antibiotic use table
    #number of days of antibiotic is randomly drawn from a truncated normal distribution
    abx_days <- round(rtnorm(ncol(los.array), mean=meanDur*timestep, sd=sdDur*timestep, lower=0))
    # Unit test - check distribution of abx distribution
    # hist(abx_days, breaks=20)
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
        abx_person <- abx_days[i]
        # Initial treatment value derived from probability, p
        rand <- runif(1, 0, 1)
        if (rand < p){
            if(abx_person > max_days){
                # if the number of generated abx is longer than the los of that person
                # have them take abx everyday for their stay
                abx.matrix[idx_end:(idx_end+max_days-1)] <- rep(1, max_days)
            }else{
                # else take abx for abx days and pad to fit max_days
                abx.matrix[idx_end:(idx_end+max_days-1)] <- c(rep(1, abx_person), rep(0, max_days-abx_person))
            }
        }else{
            # no abx taken for that person
            abx.matrix[idx_end:(idx_end+max_days-1)] <- rep(0, max_days)
        }
        # move starting position to end of previous patient
        idx_end <- idx_end+max_days
    }
    
    return(abx.matrix)
}

colo.table <- function(patient.matrix, los.array, prob_StartBact_R, prop_S_nonR){
    
    # define probabilities of importing Sensitive(S) or Resistant(R) bacteria, or low levels of sensitive (ss)
    prob_start_S <- prop_S_nonR*(1-prob_StartBact_R)
    prob_StartBact <- c(prob_start_S,prob_StartBact_R)
    
    #Generating a vector of random status with runif (change for other distribution)
    number_of_patients <- dim(los.array)[2]
    Patient_unif <- runif(number_of_patients,0,1)
    Patient_StartBact <- rep(NA, number_of_patients)
    Patient_StartBact[Patient_unif > (prob_start_S+prob_StartBact_R)] <- 'ss'
    Patient_StartBact[(Patient_unif <= (prob_start_S+prob_StartBact_R)) & (Patient_unif > prob_StartBact_R)] <- 'S'
    Patient_StartBact[Patient_unif <= prob_StartBact_R] <- 'R'
    
    #Creating array for carriage status
    array_StartBact <- matrix(NA, nrow=nrow(patient.matrix), ncol=ncol(patient.matrix))
    
    # Fill generated bacterial in the first day of each patient entering the ward
    end_idx <- 1
    for(i in 1:number_of_patients){
        array_StartBact[end_idx:(end_idx + los.array[2, i] - 1)] <- c(Patient_StartBact[i], rep(NA, los.array[2, i]-1))
        end_idx = end_idx + los.array[2, i]
    }
    
    return(array_StartBact)
}

####################4. Update values for every day  #####################
nextDay <- function(patient.matrix, los.array, abx.matrix, colo.matrix, 
                    bif, pi_ssr, repop.s1, mu_r, abx.clear, timestep=1){

    
    # adjust probabilities based on timestep
    pi_ssr <- pi_ssr/timestep
    repop.s1 <- repop.s1/timestep
    mu_r <- mu_r/timestep
    abx.clear <- abx.clear/timestep
    
    # pi_sr= probability of R transmitting to S (a proportion of pi_r if being colonised with S protects colonisation by R)
    pi_Sr <- pi_ssr - (bif*pi_ssr)
    
    # For each day (first day should be filled)
    for(i in 2:nrow(patient.matrix)){
        # Get the previous row subset of the entire matrix which represents previous day or timestep
        prev_step <- colo.matrix[i-1, ]
        # Get the indices which are already filled by a patient entering the ward for exclusion
        already_filled <- which(!is.na(colo.matrix[i, ]))
        
        # Update R
        # Get the column indices which contain R in the previous day
        R <- which(prev_step == "R")
        # Remove column indices that already have a starting bacterial state filled in
        R <- R[!(R %in% already_filled)]
        # count if there are any R on the previous day
        r_num <- length(R)
        # if there is any R (number of R > 0) in the previous day
        if(r_num){
            # Roll a random number for each R on the previous day
            roll <- runif(r_num, 0, 1)
            # All column indices which roll < mu_r (decolonization parameter) are saved to fill in as S
            decolo_idx <- R[roll < mu_r]
            # All the remaining column indices are saved to fill in as staying R the next day
            same_idx <- R[roll >= mu_r]
            # Fill in saved column as S
            colo.matrix[i, decolo_idx] <- "S"
            # Fill in saved column as R
            colo.matrix[i, same_idx] <- "R"
        }
        
        # Update S
        # Get the column indices which contain S in the previous day
        S <- which(prev_step == "S")
        # Remove column indices that already have a starting bacterial state filled in
        S <- S[!(S %in% already_filled)]
        # count if there are any S on the previous day
        s_num <- length(S)
        # if there is any S (number of S > 0) in the previous day
        if(s_num){
            # roll for transmission of R
            prob_r <- 1-((1-pi_Sr)^r_num)
            # Roll a random number for each R on the previous day for clearance
            roll_clear <- runif(s_num, 0, 1)
            # All column indices which roll < probability of clearance AND there is antibiotic used on that patient-timestep
            clear_idx <- S[abx.matrix[i-1, S] & (roll_clear < abx.clear)]
            # Clear those that pass roll and use abx to ss
            colo.matrix[i, clear_idx] <- "ss"
            # Removed those that have been cleared by abx from list of S indices
            # S <- S[!(S %in% clear_idx)] 
            # Note: Removal is not necessary since these events do not have priority. If clearance happens
            # and transmission happens, then it becomes R.
            
            # Roll a random number for each remaining S for chance of selection to
            roll_trans <- runif(length(S), 0, 1)
            r_idx <- S[roll_trans < prob_r]
            same_idx <- S[roll_trans >= prob_r]
            colo.matrix[i, r_idx] <- "R"
            colo.matrix[i, same_idx] <- "S"
        }
        
        # Update ss
        ss <- which(prev_step == "ss")
        ss <- ss[!(ss %in% already_filled)]
        if(length(ss)){
            # roll for transmission of R
            prob_r <- 1-((1-pi_ssr)^r_num)
            # roll for repop of S
            prob_s <- repop.s1
            
            # as a Gillespie approximation probability of r and s transmission should be small enough that they do not add to 1
            stopifnot((prob_r+prob_s) <= 1)
            
            roll <- runif(length(ss), 0, 1)
            r_idx <- ss[roll < prob_r]
            s_idx <- ss[(roll >= prob_r) & (roll < (prob_s+prob_r)) & !abx.matrix[i-1, ss]]
            same_idx <- ss[!(ss %in% c(r_idx, s_idx))]
            
            colo.matrix[i, s_idx] <- "S"
            colo.matrix[i, r_idx] <- "R"
            colo.matrix[i, same_idx] <- "ss"
        }
    }
    
    return(colo.matrix)
}

diff_prevalence <- function(n.bed, mean.max.los, 
                            prob_StartBact_R, prop_S_nonR, 
                            bif, pi_ssr, repop.s1, mu_r, abx.clear,
                            p, short_dur, long_dur, sdDur){
    
    old <- Sys.time() # get start time
    # DEBUG
    print(paste(n.bed, mean.max.los, 
                prob_StartBact_R, prop_S_nonR, 
                bif, pi_ssr, repop.s1, mu_r, abx.clear,
                p, short_dur, long_dur, sdDur))
    
    timestep <- 1
    n.day <- 500
    iterations <- 100
    iter_totalR <- matrix(NA, nrow = n.day, ncol = iterations)
    
    for(iter in 1:iterations){
        
        patient.matrix <- patient.table(n.bed, n.day, mean.max.los, timestep)
        los.array <- summary.los(patient.matrix)
        abx.matrix <- abx.table(patient.matrix, los.array, p, meanDur=short_dur, sdDur, timestep)
        colo.matrix <- colo.table(patient.matrix=patient.matrix, los=los.array, 
                                     prob_StartBact_R=prob_StartBact_R,prop_S_nonR=prop_S_nonR)
        
        colo_table_filled_iter <- nextDay(patient.matrix=patient.matrix, los.array=los.array, 
                                          abx.matrix=abx.matrix, colo.matrix=colo.matrix, 
                                          bif, pi_ssr, repop.s1, mu_r, abx.clear)
        
        #Summary
        df <- data.frame(colo_table_filled_iter)
        iter_totalR[, iter] <- rowSums(df == "R")    
    }
    # Discard first 1/3 runs as burn-in
    totalR_short <- mean(rowSums(iter_totalR[ceiling(n.day*1/3):nrow(iter_totalR),])/iterations/n.bed)
    
    iter_totalR <- matrix(NA, nrow = n.day, ncol = iterations)
    for(iter in 1:iterations){
        
        patient.matrix <- patient.table(n.bed, n.day, mean.max.los, timestep=1)
        los.array <- summary.los(patient.matrix)
        abx.matrix <- abx.table(patient.matrix=patient.matrix, los.array=los.array, p, meanDur=long_dur, sdDur=sdDur, timestep)
        colo.matrix <- colo.table(patient.matrix=patient.matrix, los=los.array, 
                                        prob_StartBact_R=prob_StartBact_R,prop_S_nonR=prop_S_nonR)
        
        colo_table_filled_iter <- nextDay(patient.matrix=patient.matrix, los.array=los.array, 
                                          abx.matrix=abx.matrix, colo.matrix=colo.matrix, 
                                          bif, pi_ssr, repop.s1, mu_r, abx.clear)
        
        #Summary
        df <- data.frame(colo_table_filled_iter)
        iter_totalR[, iter] <- rowSums(df == "R")    
    }
    # Discard first 1/3 runs as burn-in
    totalR_long <- mean(rowSums(iter_totalR[ceiling(n.day*1/3):nrow(iter_totalR),])/iterations/n.bed)
    
    #print(paste("totalR_long", totalR_long, "totalR_short", totalR_short))
    # print elapsed time
    new <- Sys.time() - old # calculate difference
    print(new) # print in nice format
    
    return(totalR_long - totalR_short)
}

parameters_simple<- c("n.bed", "mean.max.los", "prob_StartBact_R", "prop_S_nonR", 
               "bif", "pi_ssr", "repop.s1", "mu_r", "abx.clear", "p", "short_dur", "long_dur", "sdDur")

