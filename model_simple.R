source('msm_util_rtnorm.R')
source('los_abx_matrix.R')

#table of initial states of each patient 
colo.table <- function(patient.matrix, los.array, prob_StartBact_R, prop_S_nonR){
    
    # define probabilities of importing Sensitive(S) or Resistant(R) bacteria, or low levels of sensitive (ss)
    prob_start_S = prop_S_nonR*(1-prob_StartBact_R)
    prob_StartBact = c(prob_start_S,prob_StartBact_R)
    
    #Generating a vector of random status with runif (change for other distribution)
    number_of_patients = dim(los.array)[2]
    Patient_unif = runif(number_of_patients,0,1)
    Patient_StartBact = rep(NA, number_of_patients)
    Patient_StartBact[Patient_unif > (prob_start_S+prob_StartBact_R)] = 'ss'
    Patient_StartBact[(Patient_unif <= (prob_start_S+prob_StartBact_R)) & (Patient_unif > prob_StartBact_R)] = 'S'
    Patient_StartBact[Patient_unif <= prob_StartBact_R] = 'R'
    
    #Creating array for carriage status
    array_StartBact = matrix(NA, nrow=nrow(patient.matrix), ncol=ncol(patient.matrix))
    
    # Fill generated bacterial in the first day of each patient entering the ward
    end_idx = 1
    for(i in 1:number_of_patients){
        array_StartBact[end_idx:(end_idx + los.array[2, i] - 1)] = c(Patient_StartBact[i], rep(NA, los.array[2, i]-1))
        end_idx = end_idx + los.array[2, i]
    }
    
    return(array_StartBact)
}

####################4. Update values for every day  #####################
nextDay <- function(patient.matrix, los.array, abx.matrix, colo.matrix, 
                    bif, pi_ssr, abx.s, abx.r, repop.s1, mu_r, timestep){
    
    # adjust probabilities based on timestep
    pi_ssr = 1-(1-pi_ssr)^(1/timestep)
    repop.s1 = 1-(1-repop.s1)^(1/timestep)
    mu_r = 1-(1-mu_r)^(1/timestep)
    abx.s= 1-(1-abx.s)^(1/timestep)
    abx.r= 1-(1-abx.r)^(1/timestep)
    abx.r1 =abx.r 
    abx.r2 =abx.r 
    
    if (abx.r < 0.01) { 
        abx.matrix[abx.matrix==2]=1
    }
    
    #print(paste("pi_ssr, timestep:", pi_ssr, timestep))
    
    # pi_sr= probability of R transmitting to S (a proportion of pi_r if being colonised with S protects colonisation by R)
    pi_Sr = pi_ssr - (bif*pi_ssr)
    
    # For each day (first day should be filled)
    for(i in 2:nrow(patient.matrix)){
        # Get the previous row subset of the entire patient matrix which represents previous day or timestep
        prev_step = colo.matrix[i-1, ]
        # Get the indices which are already filled by a patient entering the ward for exclusion
        already_filled = which(!is.na(colo.matrix[i, ]))
        # Get the previous row subset of the entire abx matrix which represents previous day or timestep
        prev_abx = abx.matrix[i-1, ]
        # count if there are any R on the previous day
        R.previousday=which(prev_step == "R")
        r_num = length(R.previousday)
        
        # First scenario: those with no antibiotics on previous day 
        ##########################################
        id_noabx= which(prev_abx==0)
        
        ### Update R
        # Get the column indices which contain R in the previous day
        R = id_noabx[id_noabx %in% which(prev_step == "R")]
        # Remove column indices that already have a starting bacterial state filled in
        R = R[!(R %in% already_filled)]
        # if there is any R (number of R > 0) in the previous day
        if(length(R)){
            # Roll a random number for each R on the previous day
            roll = runif(length(R), 0, 1)
            # All column indices which roll < mu_r (decolonization parameter) are saved to fill in as S
            S_idx = R[roll < mu_r]
            # All the remaining column indices are saved to fill in as staying R the next day
            same_idx = R[roll >= mu_r]
            # Fill in saved column as S
            colo.matrix[i, S_idx] = "S"
            # Fill in saved column as R
            colo.matrix[i, same_idx] = "R"
        }
        
        ### Update S
        # Get the column indices which contain S in the previous day
        S = id_noabx[id_noabx %in% which(prev_step == "S")]
        # Remove column indices that already have a starting bacterial state filled in
        S = S[!(S %in% already_filled)]
        # if there is any S (number of S > 0) in the previous day
        if(length(S)){
            # probability for transmission of R
            prob_r = 1-((1-pi_Sr)^r_num)
            
            # Roll a random number for each S on the previous day for clearance
            roll= runif(length(S), 0, 1)
            
            r_idx = S[roll < prob_r]
            same_idx= S[roll >= prob_r]
            colo.matrix[i, r_idx] = "R"
            colo.matrix[i, same_idx] = "S"
        }
        
        ### Update ss
        ss = id_noabx[id_noabx %in% which(prev_step == "ss")]
        ss = ss[!(ss %in% already_filled)]
        if(length(ss)){
            # roll for transmission of R
            prob_r = 1-((1-pi_ssr)^r_num)
            
            roll = runif(length(ss), 0, 1)
            
            r_idx = ss[roll < prob_r]
            s_idx = ss[(roll >= prob_r) & (roll < (repop.s1+prob_r))]
            same_idx = ss[!(ss %in% c(r_idx, s_idx))]
            
            if((repop.s1+prob_r > 1)){
                stop(paste("Error stopifnot: repop.s1+prob_r >1 in First scenario: those with no antibiotics on previous day"))
            }
            
            colo.matrix[i, s_idx] = "S"
            colo.matrix[i, r_idx] = "R"
            colo.matrix[i, same_idx] = "ss"
        }
        
        # Second scenario: those with narrow antibiotics on previous day 
        ##########################################
        id_narrowabx= which(prev_abx==1)
        
        ### Update R
        # Get the column indices which contain R in the previous day
        R = id_narrowabx[id_narrowabx %in% which(prev_step == "R")]
        # Remove column indices that already have a starting bacterial state filled in
        R = R[!(R %in% already_filled)]
        # if there is any R (number of R > 0) in the previous day
        if(length(R)){
            # Fill in as R
            colo.matrix[i, R] = "R"
        }
        
        ### Update S
        # Get the column indices which contain S in the previous day
        S = id_narrowabx[id_narrowabx %in% which(prev_step == "S")]
        # Remove column indices that already have a starting bacterial state filled in
        S = S[!(S %in% already_filled)]
        # if there is any S (number of S > 0) in the previous day
        if(length(S)){
            # probability for transmission of R
            prob_r = 1-((1-pi_Sr)^r_num)
            
            # Roll a random number for each S on the previous day for clearance
            roll= runif(length(S), 0, 1)
            
            r_idx = S[roll < prob_r]
            ss_idx = S [(roll >= prob_r) & (roll < (abx.s+prob_r))]
            same_idx= S[!(S %in% c(r_idx, ss_idx))]
            
            if((abx.s+prob_r > 1)){
                stop(paste("Error stopifnot: abx.s+prob_r >1 in Second scenario: those with narrow antibiotics on previous day"))
            }
            
            colo.matrix[i, ss_idx] = "ss"
            colo.matrix[i, r_idx] = "R"
            colo.matrix[i, same_idx] = "S"
        }
        
        ### Update ss
        ss = id_narrowabx[id_narrowabx %in%which(prev_step == "ss")]
        ss = ss[!(ss %in% already_filled)]
        if(length(ss)){
            # roll for transmission of R
            prob_r = 1-((1-pi_ssr)^r_num)
            
            roll = runif(length(ss), 0, 1)
            
            r_idx = ss[roll < prob_r]
            same_idx = ss[roll >= prob_r]
            
            colo.matrix[i, r_idx] = "R"
            colo.matrix[i, same_idx] = "ss"
        }
        
        # Third scenario: those with broad antibiotics on previous day 
        ##########################################
        id_broadabx= which(prev_abx==2)
        
        ### Update R
        # Get the column indices which contain R in the previous day
        R = id_broadabx[id_broadabx %in% which(prev_step == "R")]
        # Remove column indices that already have a starting bacterial state filled in
        R = R[!(R %in% already_filled)]
        # if there is any R (number of R > 0) in the previous day
        if(length(R)){
            # Roll a random number for each R on the previous day for clearance
            roll= runif(length(R), 0, 1)
            
            # Fill in today 
            ss_idx = R[roll < abx.r2]
            same_idx = R[roll >= abx.r2]
            colo.matrix[i, ss_idx] = "ss"
            colo.matrix[i, same_idx] = "R"
        }
        
        ### Update S
        # Get the column indices which contain S in the previous day
        S = id_broadabx[id_broadabx %in% which(prev_step == "S")]
        # Remove column indices that already have a starting bacterial state filled in
        S = S[!(S %in% already_filled)]
        # if there is any S (number of S > 0) in the previous day
        if(length(S)){
            # probability for transmission of R
            prob_r = 1-((1-pi_Sr)^r_num)
            
            # Roll a random number for each S on the previous day for clearance
            roll= runif(length(S), 0, 1)
            
            r_idx = S[roll < prob_r]
            ss_idx = S [(roll >= prob_r) & (roll < (abx.r1+prob_r))]
            same_idx= S[!(S %in% c(r_idx, ss_idx))]
            
            if((abx.r1+prob_r > 1)){
                stop(paste("Error stopifnot: abx.r1+prob_r >1 in Third scenario: those with broad antibiotics on previous day "))
            }
            
            colo.matrix[i, ss_idx] = "ss"
            colo.matrix[i, r_idx] = "R"
            colo.matrix[i, same_idx] = "S"
        }
        
        ### Update ss
        ss = id_broadabx[id_broadabx %in% which(prev_step == "ss")]
        ss = ss[!(ss %in% already_filled)]
        if(length(ss)){
            # roll for transmission of R
            prob_r = 1-((1-pi_ssr)^r_num)
            
            roll = runif(length(ss), 0, 1)
            
            r_idx = ss[roll < prob_r]
            same_idx = ss[roll >= prob_r]
            
            colo.matrix[i, r_idx] = "R"
            colo.matrix[i, same_idx] = "ss"
        }
    }
    
    return(colo.matrix)
}

diff_prevalence <- function(n.bed, mean.max.los, 
                            prob_StartBact_R, prop_S_nonR, 
                            bif, pi_ssr, repop.s1, mu_r, abx.s, abx.r,
                            p.infect, cum.r.1, p.r.day1, short_dur, long_dur){
    
    old = Sys.time() # get start time
    # DEBUG
    print(paste(n.bed, mean.max.los, 
                prob_StartBact_R, prop_S_nonR, 
                bif, pi_ssr, repop.s1, mu_r, abx.s, abx.r,
                p.infect, cum.r.1, p.r.day1, short_dur, long_dur))
    
    timestep = 10
    iterations = 50
    n.day=350
    sdDur=1
    
    iter_totalR = matrix(NA, nrow = n.day*timestep, ncol = iterations)
    
    for(iter in 1:iterations){
        
        matrixes = los.abx.table(n.bed=n.bed, n.day=n.day, mean.max.los=mean.max.los, 
                                 p.infect=p.infect, p.r.day1=p.r.day1, cum.r.1=cum.r.1, 
                                 meanDur= short_dur, timestep=timestep)
        patient.matrix=matrixes[[1]]
        abx.matrix=matrixes[[2]]
        los.array = summary.los(patient.matrix=patient.matrix)
        colo.matrix = colo.table(patient.matrix=patient.matrix, los=los.array, 
                                 prob_StartBact_R=prob_StartBact_R,prop_S_nonR=prop_S_nonR)
        
        colo_table_filled_iter = nextDay(patient.matrix=patient.matrix, los.array=los.array, 
                                         abx.matrix=abx.matrix, colo.matrix=colo.matrix, 
                                         bif=bif, pi_ssr=pi_ssr, repop.s1=repop.s1, mu_r=mu_r, abx.s=abx.s, abx.r=abx.r,timestep=timestep)
        
        #Summary
        df = data.frame(colo_table_filled_iter)
        iter_totalR[, iter] = rowSums(df == "R")    
    }
    # Discard first 1/7 runs as burn-in
    totalR_short = mean(rowSums(iter_totalR[ceiling(n.day*1/7):nrow(iter_totalR), ,drop=FALSE])/iterations/n.bed)
    
    iter_totalR = matrix(NA, nrow = n.day*timestep, ncol = iterations)
    
    for(iter in 1:iterations){
        
        matrixes = los.abx.table(n.bed=n.bed, n.day=n.day, mean.max.los=mean.max.los, 
                                 p.infect=p.infect, p.r.day1=p.r.day1, cum.r.1=cum.r.1, 
                                 meanDur= long_dur, timestep=timestep)
        patient.matrix=matrixes[[1]]
        abx.matrix=matrixes[[2]]
        los.array = summary.los(patient.matrix=patient.matrix)
        colo.matrix = colo.table(patient.matrix=patient.matrix, los=los.array, 
                                 prob_StartBact_R=prob_StartBact_R,prop_S_nonR=prop_S_nonR)
        
        colo_table_filled_iter = nextDay(patient.matrix=patient.matrix, los.array=los.array, 
                                         abx.matrix=abx.matrix, colo.matrix=colo.matrix, 
                                         bif=bif, pi_ssr=pi_ssr, repop.s1=repop.s1, mu_r=mu_r, abx.s=abx.s, abx.r=abx.r,timestep=timestep)
        
        #Summary
        df = data.frame(colo_table_filled_iter)
        iter_totalR[, iter] = rowSums(df == "R")    
    }
    # Discard first 1/7 runs as burn-in
    totalR_long = mean(rowSums(iter_totalR[ceiling(n.day*1/7):nrow(iter_totalR), ,drop=FALSE])/iterations/n.bed)
    
    #print(paste("totalR_long", totalR_long, "totalR_short", totalR_short))
    # print elapsed time
    new = Sys.time() - old # calculate difference
    print(new) # print in nice format
    
    return(totalR_long - totalR_short)
}

parameters_simple<- c("n.bed", "mean.max.los", 
                      "prob_StartBact_R", "prop_S_nonR", 
                      "bif", "pi_ssr", "repop.s1", "mu_r", 
                      "abx.s", "abx.r", "p.infect", "cum.r.1", 'p.r.day1', "short_dur", "long_dur")
