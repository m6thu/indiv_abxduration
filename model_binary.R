source('msm_util_rtnorm.R')
source('los_abx_matrix.R')

#(1 for short duration and 1 for long duration) 
#simulate inpatients with various lengths of stay
#allocate various duration of antibiotics for each patient

#################Generate baseline carriage status 

colo.table <- function(patient.matrix, los.array, prob_StartBact_R, prop_S_nonR, prop_Sr_inR, prop_sr_inR){
    
    prob_start_S = prop_S_nonR*(1-prob_StartBact_R)
    prob_start_ss = 1-prob_start_S-prob_StartBact_R
    prob_start_Sr = prop_Sr_inR*prob_StartBact_R
    prob_start_sr = prop_sr_inR*prob_StartBact_R
    prob_start_sR = prob_StartBact_R-prob_start_Sr-prob_start_sr
    prob_StartBact_bi = c(prob_start_S,prob_start_Sr,prob_start_sR,prob_start_sr)

    #Generating a vector of random status with runif (change for other distribution)
    number_of_patients = dim(los.array)[2]
    Patient_unif = runif(number_of_patients,0,1)
    Patient_StartBact = rep(NA, number_of_patients)
    Patient_StartBact[Patient_unif > sum(prob_StartBact_bi)] = 'ss'
    Patient_StartBact[(Patient_unif > sum(prob_StartBact_bi[1:3])) & (Patient_unif <= sum(prob_StartBact_bi))] = 'sr'
    Patient_StartBact[(Patient_unif > sum(prob_StartBact_bi[1:2])) & (Patient_unif <= sum(prob_StartBact_bi[1:3]))] = 'sR'
    Patient_StartBact[(Patient_unif > sum(prob_StartBact_bi[1])) & (Patient_unif <= sum(prob_StartBact_bi[1:2]))] = 'Sr'
    Patient_StartBact[Patient_unif <= prob_StartBact_bi[1]] = 'S'
    
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

#################### Update values for every day  
nextDay <- function(patient.matrix, abx.matrix, colo.matrix, 
                        pi_ssr, bif, mu1, mu2, mu_r, repop.r1, repop.r2,
                        repop.s1, repop.s2, abx.r, abx.s, timestep){
    
    # adjust probabilities based on timestep
    pi_ssr = 1-(1-pi_ssr)^(1/timestep)
    mu1 = 1-(1-mu1)^(1/timestep)
    mu2 = 1-(1-mu2)^(1/timestep)
    mu_r=1-(1-mu_r)^(1/timestep)
    repop.r1 = 1-(1-repop.r1)^(1/timestep)
    repop.r2 = 1-(1-repop.r2)^(1/timestep)
    repop.s1 = 1-(1-repop.s1)^(1/timestep)
    repop.s2 = 1-(1-repop.s2)^(1/timestep)
    abx.r = 1-(1-abx.r)^(1/timestep)
    abx.s = 1-(1-abx.s)^(1/timestep)
    abx.r1 =abx.r 
    abx.r2 =abx.r 
    
    if (abx.r < 0.01) { 
        abx.matrix[abx.matrix==2]=1
    }
    
    pi_r1 = pi_ssr- (pi_ssr * bif)                 # pi_ssr= probability of R transmitting to s to become sr 
    
    # For each day
    for(i in 2:nrow(patient.matrix)){ # start from 2 because for each day (first day should be filled)
        # Get the previous row subset of the entire matrix which represents previous day or timestep
        prev_step = colo.matrix[i-1, ]
        # Get the indices which are already filled by a patient entering the ward for exclusion
        already_filled = which(!is.na(colo.matrix[i, ]))
        # Get the previous row subset of the entire abx matrix which represents previous day or timestep
        prev_abx = abx.matrix[i-1, ]
        # count if there are any sR on the previous day
        R.previousday=which(prev_step == "sR")
        r_num = length(R.previousday)
        
        
        # First scenario: those with no antibiotics on previous day 
        ##########################################
        id_noabx= which(prev_abx==0)
        
        # Update S
        # Get the column indices which contain S in the previous day
        S = id_noabx[id_noabx %in% which(prev_step == "S")]
        # Remove column indices that already have a starting bacterial state filled in
        S = S[!(S %in% already_filled)]
        # if there is any S (number of S > 0) in the previous day
        if(length(S)){
            # roll for transmission of R
            prob_r = 1-((1-pi_r1)^r_num)
            
            # Roll a random number for each S on the previous day for clearance
            roll= runif(length(S), 0, 1)
            
            Sr_idx = S[roll < prob_r]
            same_idx= S[roll >= prob_r]
            colo.matrix[i, Sr_idx] = "Sr"
            colo.matrix[i, same_idx] = "S"
        }
        
        # Update ss
        ss = id_noabx[id_noabx %in% which(prev_step == "ss")]
        ss = ss[!(ss %in% already_filled)]
        if(length(ss)){
            # roll for transmission of R
            prob_r = 1-((1-pi_ssr)^r_num)
            
            roll = runif(length(ss), 0, 1)
            
            ssr_idx = ss[roll < prob_r]
            s_idx = ss[(roll >= prob_r) & (roll < (repop.s1+prob_r))]
            same_idx = ss[!(ss %in% c(ssr_idx, s_idx))]
            
            if((repop.s1+prob_r > 1)){
                stop(paste("Error stopifnot: repop.s1+prob_r >1 in First scenario: those with no antibiotics on previous day"))
            }
            
            colo.matrix[i, s_idx] = "S"
            colo.matrix[i, ssr_idx] = "sR"
            colo.matrix[i, same_idx] = "ss"
        }
        
        # Update Sr
        Sr = id_noabx[id_noabx %in% which(prev_step == "Sr")]
        Sr = Sr[!(Sr %in% already_filled)]
        
        if(length(Sr)){
            
            # Roll a random number for each S on the previous day for clearance
            roll= runif(length(Sr), 0, 1)
            
            s_idx = Sr[roll < mu1]
            same_idx= Sr[roll >= mu1]
            colo.matrix[i, s_idx] = "S"
            colo.matrix[i, same_idx] = "Sr"
        }
        
        # Update sr
        sr = id_noabx[id_noabx %in% which(prev_step == "sr")]
        sr = sr[!(sr %in% already_filled)]
        
        if(length(sr)){
            
            if(repop.r1+repop.s2 + mu2 > 1){
                stop(paste("Error stopifnot: repop.r1+repop.s2 + mu2 > 1 in First scenario: those with no antibiotics on previous day"))
            }
            
            roll= runif(length(sr), 0, 1)
            
            # If roll passes repop.r1 R grows
            sR_idx = sr[roll < repop.r1]
            # if roll does not pass repop.r1 event and passes repop.s2 
            Sr_idx = sr[roll > repop.r1  & roll <= repop.s2+repop.r1]
            # if roll does not pass repop.r1 and repop.s2, and passes mu2
            ss_idx = sr[roll>repop.s2+repop.r1 & roll <= repop.s2+repop.r1+mu2]
            
            same_idx= sr[!(sr %in% c(sR_idx, Sr_idx, ss_idx))]
            
            colo.matrix[i, sR_idx] = "sR"
            colo.matrix[i, Sr_idx] = "Sr"
            colo.matrix[i, ss_idx] = "ss"
            colo.matrix[i, same_idx] = "sr"
        }
        
        # Update sR
        sR = id_noabx[id_noabx %in% which(prev_step == "sR")]
        sR = sR[!(sR %in% already_filled)]
        if(length(sR)){
            
            roll = runif(length(sR), 0, 1)
            
            ssr_idx = sR[roll < mu_r]
            same_idx = sR[roll >= mu_r]
            
            colo.matrix[i, ssr_idx] = "sr"
            colo.matrix[i, same_idx] = "sR"
        }
        
        # Second scenario: those with antibiotics.s on previous day 
        ##########################################
        id_narrowabx= which(prev_abx==1)
        
        # Update S
        # Get the column indices which contain S in the previous day
        S = id_narrowabx[id_narrowabx %in% which(prev_step == "S")]
        # Remove column indices that already have a starting bacterial state filled in
        S = S[!(S %in% already_filled)]
        # if there is any S (number of S > 0) in the previous day
        if(length(S)){
            
            # probability of transmission 
            prob_r = 1-((1-pi_r1)^r_num)
            
            if((abx.s+prob_r > 1)){
                stop(paste("Error stopifnot: abx.s+prob_r >1 in Second scenario: those with narrow antibiotics on previous day"))
            }
            
            # Roll a random number for each S on the previous day for clearance
            roll= runif(length(S), 0, 1)
            
            Sr_idx = S[roll < prob_r]
            ss_idx = S[(roll >= prob_r) & (roll < (abx.s+prob_r))]
            same_idx= S[!(S %in% c(Sr_idx, ss_idx))]
            
            colo.matrix[i, ss_idx] = "ss"
            colo.matrix[i, Sr_idx] = "Sr"
            colo.matrix[i, same_idx] = "S"
        }
        
        # Update ss
        ss = id_narrowabx[id_narrowabx %in% which(prev_step == "ss")]
        ss = ss[!(ss %in% already_filled)]
        if(length(ss)){
            # roll for transmission of R
            prob_r = 1-((1-pi_ssr)^r_num)
            
            roll = runif(length(ss), 0, 1)
            
            ssr_idx = ss[roll < prob_r]
            same_idx = ss[roll >= prob_r]
            
            colo.matrix[i, ssr_idx] = "sR"
            colo.matrix[i, same_idx] = "ss"
        }
        
        # Update Sr
        Sr = id_narrowabx[id_narrowabx %in% which(prev_step == "Sr")]
        Sr = Sr[!(Sr %in% already_filled)]
        
        if(length(Sr)){
            
            # Roll a random number for each Sr 
            roll= runif(length(Sr), 0, 1)
            
            sr_idx = Sr[roll < abx.s]
            same_idx= Sr[roll >= abx.s]
            
            colo.matrix[i, sr_idx] = "sr"
            colo.matrix[i, same_idx] = "Sr"
        }
        
        # Update sr
        sr = id_narrowabx[id_narrowabx %in% which(prev_step == "sr")]
        sr = sr[!(sr %in% already_filled)]
        
        if(length(sr)){
            
            roll= runif(length(sr), 0, 1)
            
            sR_idx = sr[roll < repop.r2]
            same_idx= sr[roll >= repop.r2]
            
            colo.matrix[i, sR_idx] = "sR"
            colo.matrix[i, same_idx] = "sr"
        }
        
        # Update sR 
        sR = id_narrowabx[id_narrowabx %in% which(prev_step == "sR")]
        sR = sR[!(sR %in% already_filled)]
        
        if(length(sR)){
            
            colo.matrix[i, sR] = "sR"
        }
        
        # Third scenario: those with antibiotics.r on previous day 
        ##########################################
        id_broadabx= which(prev_abx==2)
        
        # Update S
        # Get the column indices which contain S in the previous day
        S = id_broadabx[id_broadabx %in% which(prev_step == "S")]
        # Remove column indices that already have a starting bacterial state filled in
        S = S[!(S %in% already_filled)]
        # if there is any S (number of S > 0) in the previous day
        if(length(S)){
            # roll for transmission of R
            prob_r = 1-((1-pi_r1)^r_num)
            
            if((abx.r1+prob_r > 1)){
                stop(paste("Error stopifnot: abx.r1+prob_r >1 in Third scenario: those with broad antibiotics on previous day "))
            }
            
            # Roll a random number for each S on the previous day for clearance
            roll= runif(length(S), 0, 1)
            
            Sr_idx = S[roll < prob_r]
            ss_idx = S[(roll >= prob_r) & (roll < (abx.r1+prob_r))]
            same_idx= S[!(S %in% c(Sr_idx, ss_idx))]
            
            colo.matrix[i, ss_idx] = "ss"
            colo.matrix[i, Sr_idx] = "Sr"
            colo.matrix[i, same_idx] = "S"
        }
        
        # Update ss
        ss = id_broadabx[id_broadabx %in% which(prev_step == "ss")]
        ss = ss[!(ss %in% already_filled)]
        if(length(ss)){
            # roll for transmission of R
            prob_r = 1-((1-pi_ssr)^r_num)
            
            roll = runif(length(ss), 0, 1)
            
            ssr_idx = ss[roll < prob_r]
            same_idx = ss[roll >= prob_r]
            
            colo.matrix[i, ssr_idx] = "sR"
            colo.matrix[i, same_idx] = "ss"
        }
        
        # Update Sr
        Sr = id_broadabx[id_broadabx %in% which(prev_step == "Sr")]
        Sr = Sr[!(Sr %in% already_filled)]
        
        if(length(Sr)){
            
            # Roll a random number for each S on the previous day for clearance
            roll= runif(length(Sr), 0, 1)
            
            sr_idx = Sr[roll < abx.r1]
            same_idx= Sr[roll >= abx.r1]
            colo.matrix[i, sr_idx] = "sr"
            colo.matrix[i, same_idx] = "Sr"
        }
        
        # Update sr
        sr = id_broadabx[id_broadabx %in% which(prev_step == "sr")]
        sr = sr[!(sr %in% already_filled)]
        
        if(length(sr)){
            
            roll= runif(length(sr), 0, 1)
            
            ss_idx = sr[roll < abx.r2]
            same_idx= sr[roll >= abx.r2]
            colo.matrix[i, ss_idx] = "ss"
            colo.matrix[i, same_idx] = "sr"
        }
        
        # Update sR
        sR = id_broadabx[id_broadabx %in% which(prev_step == "sR")]
        sR = sR[!(sR %in% already_filled)]
        if(length(sR)){
            
            roll = runif(length(sR), 0, 1)
            
            sr_idx = sR[roll < abx.r2 + mu_r]
            same_idx= sR[roll >= abx.r2 + mu_r]
            colo.matrix[i, sr_idx] = "sr"
            colo.matrix[i, same_idx] = "sR"
        }
    }
    
    return(colo.matrix)
}

diff_prevalence <- function(n.bed, mean.max.los, 
                            prob_StartBact_R, prop_S_nonR, prop_Sr_inR, prop_sr_inR,
                            bif, pi_ssr, repop.s1, repop.s2, repop.r1, repop.r2,
                            mu1, mu2, mu_r, abx.s, abx.r, 
                            p.infect, cum.r.1, p.r.day1, short_dur, long_dur){
    
    old = Sys.time() # get start time
    # DEBUG
    print(paste(n.bed, mean.max.los, 
                prob_StartBact_R, prop_S_nonR, prop_Sr_inR, prop_sr_inR,
                bif, pi_ssr, repop.s1, repop.s2, repop.r1, repop.r2,
                mu1, mu2, mu_r, abx.s, abx.r, 
                p.infect, cum.r.1, p.r.day1, short_dur, long_dur))
    
    timestep = 3
    n.day = 300
    iterations = 2
    
    iter_totalsR = matrix(NA, nrow = n.day, ncol = iterations)
    iter_totalr_or_R= matrix(NA, nrow = n.day, ncol = iterations)
    for(iter in 1:iterations){
        matrixes = los.abx.table(n.bed=n.bed, n.day=n.day, mean.max.los=mean.max.los, 
                                 p.infect=p.infect, p.r.day1=p.r.day1, cum.r.1=cum.r.1, 
                                 meanDur= short_dur, timestep=timestep)
        patient.matrix=matrixes[[1]]
        abx.matrix=matrixes[[2]]
        los.array = summary.los(patient.matrix=patient.matrix)
        colo.matrix = colo.table(patient.matrix=patient.matrix, los=los.array, 
                                  prob_StartBact_R=prob_StartBact_R, prop_S_nonR=prop_S_nonR, prop_Sr_inR=prop_Sr_inR, prop_sr_inR=prop_sr_inR)
        
        colo.matrix_filled_iter = nextDay(patient.matrix=patient.matrix, abx.matrix=abx.matrix, colo.matrix=colo.matrix, 
                                          pi_ssr=pi_ssr, bif=bif, mu1=mu1, mu2=mu2, mu_r=mu_r, repop.r1=repop.r1, repop.r2=repop.r2,
                                          repop.s1=repop.s1, repop.s2=repop.s2, abx.r=abx.r, abx.s=abx.s, timestep=timestep)
        
        #Summary
        df = data.frame(colo.matrix_filled_iter)
        iter_totalsR[, iter] = rowMeans(matrix(rowSums(df == "sR"), ncol=timestep, byrow=T))
        iter_totalr_or_R[, iter] = rowMeans(matrix(rowSums(df == "sR" | df == "sr" |df == "Sr"), ncol=timestep, byrow=T))
        #print("end iteration loop")
    }
    
    totalsR_short = mean(rowSums(iter_totalsR[151:nrow(iter_totalsR),, drop=FALSE])/iterations/n.bed)
    totalr_or_R_short = mean(rowSums(iter_totalr_or_R[151:nrow(iter_totalr_or_R),, drop=FALSE])/iterations/n.bed)
    
    iter_totalsR = matrix(NA, nrow = n.day, ncol = iterations)
    iter_totalr_or_R = matrix(NA, nrow = n.day, ncol = iterations)
    
    for(iter in 1:iterations){
        matrixes = los.abx.table(n.bed=n.bed, n.day=n.day, mean.max.los=mean.max.los, 
                                 p.infect=p.infect, p.r.day1=p.r.day1, cum.r.1=cum.r.1, 
                                 meanDur= long_dur, timestep=timestep)
        patient.matrix=matrixes[[1]]
        abx.matrix=matrixes[[2]]
        los.array = summary.los(patient.matrix=patient.matrix)
        colo.matrix = colo.table(patient.matrix=patient.matrix, los=los.array, 
                                  prob_StartBact_R, prop_S_nonR, prop_Sr_inR, prop_sr_inR)
        
        colo.matrix_filled_iter = nextDay(patient.matrix= patient.matrix, abx.matrix=abx.matrix, colo.matrix=colo.matrix, 
                                          pi_ssr=pi_ssr, bif=bif, mu1=mu1, mu2=mu2, mu_r=mu_r, repop.r1=repop.r1, repop.r2=repop.r2,
                                          repop.s1=repop.s1, repop.s2=repop.s2, abx.r=abx.r, abx.s=abx.s, timestep=timestep)
        #Summary
        df = data.frame(colo.matrix_filled_iter)
        iter_totalsR[,iter] = rowMeans(matrix(rowSums(df == "sR"), ncol=timestep, byrow=T))
        iter_totalr_or_R[, iter] = rowMeans(matrix(rowSums(df == "sR" | df == "sr" |df == "Sr"), ncol=timestep, byrow=T))
        #print("end iteration loop")
    }
    totalsR_long = mean(rowSums(iter_totalsR[151:nrow(iter_totalsR),, drop=FALSE])/iterations/n.bed)
    totalr_or_R_long = mean(rowSums(iter_totalr_or_R[151:nrow(iter_totalr_or_R),, drop=FALSE])/iterations/n.bed)
    
    #print(paste("totalsR_long", totalsR_long, "totalsR_short", totalsR_short))
    # print elapsed time
    new = Sys.time() - old # calculate difference
    print(new) # print in nice format
    
    return(array(c(totalsR_long, totalsR_short,totalsR_long - totalsR_short)))
}

prevalence <- function(n.bed, mean.max.los, 
                            prob_StartBact_R, prop_S_nonR, prop_Sr_inR, prop_sr_inR,
                            bif, pi_ssr, repop.s1, repop.s2, repop.r1, repop.r2,
                            mu1, mu2, mu_r, abx.s, abx.r, 
                            p.infect, cum.r.1, p.r.day1, meanDur){
    
    old = Sys.time() # get start time

    timestep = 3
    n.day = 300
    iterations = 125
    
    iter_totalsR = matrix(NA, nrow = n.day, ncol = iterations)
    for(iter in 1:iterations){
        matrixes = los.abx.table(n.bed=n.bed, n.day=n.day, mean.max.los=mean.max.los, 
                                 p.infect=p.infect, p.r.day1=p.r.day1, cum.r.1=cum.r.1, 
                                 meanDur= meanDur, timestep=timestep)
        patient.matrix=matrixes[[1]]
        abx.matrix=matrixes[[2]]
        los.array = summary.los(patient.matrix=patient.matrix)
        colo.matrix = colo.table(patient.matrix=patient.matrix, los=los.array, 
                                 prob_StartBact_R=prob_StartBact_R, prop_S_nonR=prop_S_nonR, prop_Sr_inR=prop_Sr_inR, prop_sr_inR=prop_sr_inR)
        
        colo.matrix_filled_iter = nextDay(patient.matrix=patient.matrix, abx.matrix=abx.matrix, colo.matrix=colo.matrix, 
                                          pi_ssr=pi_ssr, bif=bif, mu1=mu1, mu2=mu2, mu_r=mu_r, repop.r1=repop.r1, repop.r2=repop.r2,
                                          repop.s1=repop.s1, repop.s2=repop.s2, abx.r=abx.r, abx.s=abx.s, timestep=timestep)
        
        #Summary
        df = data.frame(colo.matrix_filled_iter)
        iter_totalsR[, iter] = rowMeans(matrix(rowSums(df == "sR"), ncol=timestep, byrow=T))
        #print("end iteration loop")
    }
    
    totalsR = mean(rowSums(iter_totalsR[151:nrow(iter_totalsR),, drop=FALSE])/iterations/n.bed)
    
    #print(paste("totalsR_long", totalsR_long, "totalsR_short", totalsR_short))
    # print elapsed time
    new = Sys.time() - old # calculate difference
    print(new) # print in nice format
    
    return(totalsR)
}

parameters_prevalence_binary <- c("n.bed", "mean.max.los", 
                                       "prob_StartBact_R", "prop_S_nonR", "prop_Sr_inR", "prop_sr_inR",
                                       "bif", "pi_ssr", "repop.s1", "repop.s2", "repop.r1","repop.r2",
                                       "mu1", "mu2", "mu_r", "abx.s", "abx.r",
                                       "p.infect", "cum.r.1", "p.r.day1",
                                       "meanDur")

parameters_diff_prevalence_binary <- c("n.bed", "mean.max.los", 
                      "prob_StartBact_R", "prop_S_nonR", "prop_Sr_inR", "prop_sr_inR",
                      "bif", "pi_ssr", "repop.s1", "repop.s2", "repop.r1","repop.r2",
                      "mu1", "mu2", "mu_r", "abx.s", "abx.r",
                      "p.infect", "cum.r.1", "p.r.day1",
                      "short_dur", "long_dur")



