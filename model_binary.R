source('msm_util_rtnorm.R')
source('los_abx_matrix.R')

#(1 for short duration and 1 for long duration) 
#simulate inpatients with various lengths of stay
#allocate various duration of antibiotics for each patient

#################Generate baseline carriage status 

colo.table <- function(patient.matrix, los.array, prop_R, prop_S_nonR, prop_Sr_inR, prop_sr_inR){
    
    prob_start_S = prop_S_nonR*(1-prop_R)
    prob_start_ss = 1-prob_start_S-prop_R
    prob_start_Sr = prop_Sr_inR*prop_R
    prob_start_sr = prop_sr_inR*prop_R
    prob_start_sR = prop_R-prob_start_Sr-prob_start_sr
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
                        pi_ssr, bif, mu, repop.r, repop.s, abx.r, abx.s, timestep){
    
    # adjust probabilities based on timestep
    pi_ssr = 1-(1-pi_ssr)^(1/timestep)
    mu = 1-(1-mu)^(1/timestep)
    repop.r = 1-(1-repop.r)^(1/timestep)
    repop.s = 1-(1-repop.s)^(1/timestep)
    abx.r = 1-(1-abx.r)^(1/timestep)
    abx.s = 1-(1-abx.s)^(1/timestep)
    
    if (abx.r < 0.05) { #if abx.r is ineffective in scenario B - resistance = CRE 
      abx.r.aginst.s =abx.s #remains effective for s  
      abx.r.aginst.r =abx.r #ineffective for r 
    } else { # if abx.r is effective in scenario A - resistance = ESBL
      abx.r.aginst.s =abx.r #effective for both s and r 
      abx.r.aginst.r =abx.r
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
            prop_R = 1-((1-pi_r1)^r_num)
            
            # Roll a random number for each S on the previous day for clearance
            roll= runif(length(S), 0, 1)
            
            Sr_idx = S[roll < prop_R]
            same_idx= S[roll >= prop_R]
            colo.matrix[i, Sr_idx] = "Sr"
            colo.matrix[i, same_idx] = "S"
        }
        
        # Update ss
        ss = id_noabx[id_noabx %in% which(prev_step == "ss")]
        ss = ss[!(ss %in% already_filled)]
        if(length(ss)){
            # roll for transmission of R
            prop_R = 1-((1-pi_ssr)^r_num)
            
            roll = runif(length(ss), 0, 1)
            
            ssr_idx = ss[roll < prop_R]
            s_idx = ss[(roll >= prop_R) & (roll < (repop.s+prop_R))]
            same_idx = ss[!(ss %in% c(ssr_idx, s_idx))]
            
            if((repop.s+prop_R > 1)){
                stop(paste("Error stopifnot: repop.s+prop_R >1 in First scenario: those with no antibiotics on previous day"))
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
            
            s_idx = Sr[roll < mu]
            same_idx= Sr[roll >= mu]
            colo.matrix[i, s_idx] = "S"
            colo.matrix[i, same_idx] = "Sr"
        }
        
        # Update sr
        sr = id_noabx[id_noabx %in% which(prev_step == "sr")]
        sr = sr[!(sr %in% already_filled)]
        
        if(length(sr)){
            
            if(repop.r+repop.s + mu > 1){
                stop(paste("Error stopifnot: repop.r+repop.s + mu > 1 in First scenario: those with no antibiotics on previous day"))
            }
            
            roll= runif(length(sr), 0, 1)
            
            # If roll passes repop.r R grows
            sR_idx = sr[roll < repop.r]
            # if roll does not pass repop.r event and passes repop.s 
            Sr_idx = sr[roll > repop.r  & roll <= repop.s+repop.r]
            # if roll does not pass repop.r and repop.s, and passes mu
            ss_idx = sr[roll>repop.s+repop.r & roll <= repop.s+repop.r+mu]
            
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
            
            ssr_idx = sR[roll < mu]
            same_idx = sR[roll >= mu]
            
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
            prop_R = 1-((1-pi_r1)^r_num)
            
            if((abx.s+prop_R > 1)){
                stop(paste("Error stopifnot: abx.s+prop_R >1 in Second scenario: those with narrow antibiotics on previous day"))
            }
            
            # Roll a random number for each S on the previous day for clearance
            roll= runif(length(S), 0, 1)
            
            Sr_idx = S[roll < prop_R]
            ss_idx = S[(roll >= prop_R) & (roll < (abx.s+prop_R))]
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
            prop_R = 1-((1-pi_ssr)^r_num)
            
            roll = runif(length(ss), 0, 1)
            
            ssr_idx = ss[roll < prop_R]
            same_idx = ss[roll >= prop_R]
            
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
            
            sR_idx = sr[roll < repop.r]
            same_idx= sr[roll >= repop.r]
            
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
            prop_R = 1-((1-pi_r1)^r_num)
            
            if((abx.r.aginst.s+prop_R > 1)){
                stop(paste("Error stopifnot: abx.r.aginst.s+prop_R >1 in Third scenario: those with broad antibiotics on previous day "))
            }
            
            # Roll a random number for each S on the previous day for clearance
            roll= runif(length(S), 0, 1)
            
            Sr_idx = S[roll < prop_R]
            ss_idx = S[(roll >= prop_R) & (roll < (abx.r.aginst.s+prop_R))]
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
            prop_R = 1-((1-pi_ssr)^r_num)
            
            roll = runif(length(ss), 0, 1)
            
            ssr_idx = ss[roll < prop_R]
            same_idx = ss[roll >= prop_R]
            
            colo.matrix[i, ssr_idx] = "sR"
            colo.matrix[i, same_idx] = "ss"
        }
        
        # Update Sr
        Sr = id_broadabx[id_broadabx %in% which(prev_step == "Sr")]
        Sr = Sr[!(Sr %in% already_filled)]
        
        if(length(Sr)){
            
            # Roll a random number for each S on the previous day for clearance
            roll= runif(length(Sr), 0, 1)
            
            sr_idx = Sr[roll < abx.r.aginst.s]
            same_idx= Sr[roll >= abx.r.aginst.s]
            
            colo.matrix[i, sr_idx] = "sr"
            colo.matrix[i, same_idx] = "Sr"
        }
        
        # Update sr
        sr = id_broadabx[id_broadabx %in% which(prev_step == "sr")]
        sr = sr[!(sr %in% already_filled)]
        
        if(length(sr)){
            
            roll= runif(length(sr), 0, 1)
            
            ss_idx = sr[roll < abx.r.aginst.r]
            same_idx= sr[roll >= abx.r.aginst.r]
            colo.matrix[i, ss_idx] = "ss"
            colo.matrix[i, same_idx] = "sr"
        }
        
        # Update sR
        sR = id_broadabx[id_broadabx %in% which(prev_step == "sR")]
        sR = sR[!(sR %in% already_filled)]
        if(length(sR)){
            
            roll = runif(length(sR), 0, 1)
            
            sr_idx = sR[roll < abx.r.aginst.r]
            same_idx= sR[roll >= abx.r.aginst.r]
            
            colo.matrix[i, sr_idx] = "sr"
            colo.matrix[i, same_idx] = "sR"
        }
    }
    
    return(colo.matrix)
}

diff_prevalence <- function(n.bed, max.los, 
                            prop_R, prop_S_nonR, prop_Sr_inR, prop_sr_inR,
                            bif, pi_ssr, repop.s, repop.r,
                            mu, abx.s, abx.r, 
                            p.infect, cum.r.1, p.r.day1, short_dur, long_dur){
    
    old = Sys.time() # get start time
    # DEBUG
    print(paste(n.bed, max.los, 
                prop_R, prop_S_nonR, prop_Sr_inR, prop_sr_inR,
                bif, pi_ssr, repop.s, repop.r,
                mu, abx.s, abx.r, 
                p.infect, cum.r.1, p.r.day1, short_dur, long_dur))
    
    timestep = 1
    n.day = 300
    
    if (abx.r>0.01){ #scenario A
      iterations= 125
    } else { #scenario B
      iterations= 100
    }
    
    iter_totalsR = matrix(NA, nrow = n.day, ncol = iterations)
    iter_totalr_or_R= matrix(NA, nrow = n.day, ncol = iterations)
    for(iter in 1:iterations){
        matrixes = los.abx.table(n.bed=n.bed, n.day=n.day, max.los=max.los, 
                                 p.infect=p.infect, p.r.day1=p.r.day1, cum.r.1=cum.r.1, 
                                 meanDur= short_dur, timestep=timestep)
        patient.matrix=matrixes[[1]]
        abx.matrix=matrixes[[2]]
        los.array = summary.los(patient.matrix=patient.matrix)
        colo.matrix = colo.table(patient.matrix=patient.matrix, los=los.array, 
                                  prop_R=prop_R, prop_S_nonR=prop_S_nonR, prop_Sr_inR=prop_Sr_inR, prop_sr_inR=prop_sr_inR)
        
        colo.matrix_filled_iter = nextDay(patient.matrix=patient.matrix, abx.matrix=abx.matrix, colo.matrix=colo.matrix, 
                                          pi_ssr=pi_ssr, bif=bif, mu=mu, repop.r=repop.r, 
                                          repop.s=repop.s, abx.r=abx.r, abx.s=abx.s, timestep=timestep)
        
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
        matrixes = los.abx.table(n.bed=n.bed, n.day=n.day, max.los=max.los, 
                                 p.infect=p.infect, p.r.day1=p.r.day1, cum.r.1=cum.r.1, 
                                 meanDur= long_dur, timestep=timestep)
        patient.matrix=matrixes[[1]]
        abx.matrix=matrixes[[2]]
        los.array = summary.los(patient.matrix=patient.matrix)
        colo.matrix = colo.table(patient.matrix=patient.matrix, los=los.array, 
                                  prop_R, prop_S_nonR, prop_Sr_inR, prop_sr_inR)
        
        colo.matrix_filled_iter = nextDay(patient.matrix= patient.matrix, abx.matrix=abx.matrix, colo.matrix=colo.matrix, 
                                          pi_ssr=pi_ssr, bif=bif, mu=mu, repop.r=repop.r,
                                          repop.s=repop.s, abx.r=abx.r, abx.s=abx.s, timestep=timestep)
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

prevalence <- function(n.bed, max.los, 
                            prop_R, prop_S_nonR, prop_Sr_inR, prop_sr_inR,
                            bif, pi_ssr, repop.s, repop.r,
                            mu, abx.s, abx.r, 
                            p.infect, cum.r.1, p.r.day1, meanDur){
    
    old = Sys.time() # get start time
    timestep = 1
    n.day = 300
    
    if (abx.r>0.01){ #scenario A
      iterations= 125
    } else { #scenario B
      iterations= 100
    }
    
    iter_totalsR = matrix(NA, nrow = n.day, ncol = iterations)
    for(iter in 1:iterations){
        matrixes = los.abx.table(n.bed=n.bed, n.day=n.day, max.los=max.los, 
                                 p.infect=p.infect, p.r.day1=p.r.day1, cum.r.1=cum.r.1, 
                                 meanDur= meanDur, timestep=timestep)
        patient.matrix=matrixes[[1]]
        abx.matrix=matrixes[[2]]
        los.array = summary.los(patient.matrix=patient.matrix)
        colo.matrix = colo.table(patient.matrix=patient.matrix, los=los.array, 
                                 prop_R=prop_R, prop_S_nonR=prop_S_nonR, prop_Sr_inR=prop_Sr_inR, prop_sr_inR=prop_sr_inR)
        
        colo.matrix_filled_iter = nextDay(patient.matrix=patient.matrix, abx.matrix=abx.matrix, colo.matrix=colo.matrix, 
                                          pi_ssr=pi_ssr, bif=bif, mu=mu, repop.r=repop.r, 
                                          repop.s=repop.s, abx.r=abx.r, abx.s=abx.s, timestep=timestep)
        
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

parameters_prevalence_binary <- c("n.bed", "max.los", 
                                       "prop_R", "prop_S_nonR", "prop_Sr_inR", "prop_sr_inR",
                                       "bif", "pi_ssr", "repop.s", "repop.r",
                                       "mu", "abx.s", "abx.r",
                                       "p.infect", "cum.r.1", "p.r.day1",
                                       "meanDur")

parameters_diff_prevalence_binary <- c("n.bed", "max.los", 
                      "prop_R", "prop_S_nonR", "prop_Sr_inR", "prop_sr_inR",
                      "bif", "pi_ssr", "repop.s","repop.r",
                      "mu", "abx.s", "abx.r",
                      "p.infect", "cum.r.1", "p.r.day1",
                      "short_dur", "long_dur")



