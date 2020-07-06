source('msm_util_rtnorm.R')
source('los_abx_matrix.R')

#table of initial states of each patient 
colo.table <- function(patient.matrix, los.array, prop_R, prop_S){
  
  # define probabilities of importing Sensitive(S) or Resistant(R) bacteria, or low levels of sensitive (ss)
  prop_start_S = prop_S*(1-prop_R)
  
  #Generating a vector of random status with runif (change for other distribution)
  number_of_patients = dim(los.array)[2]
  Patient_unif = runif(number_of_patients, 0, 1)
  Patient_StartBact = rep(NA, number_of_patients)
  Patient_StartBact[Patient_unif > prop_start_S + prop_R] = 'ss'
  Patient_StartBact[Patient_unif <= prop_start_S + prop_R & Patient_unif > prop_R] = 'S'
  Patient_StartBact[Patient_unif <= prop_R] = 'R'
  
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
                    bif, pi_ssr, abx.s, abx.r, repop.s, mu, timestep){
  
  n.bed = ncol(patient.matrix)
  
  # adjust probabilities based on timestep
  if (timestep > 1) {
    pi_ssr = 1 - (1 - pi_ssr) ^ (1 / timestep)
    repop.s = 1 - (1 - repop.s) ^ (1 / timestep)
    mu = 1 - (1 - mu) ^ (1 / timestep)
    abx.s = 1 - (1 - abx.s) ^ (1 / timestep)
    abx.r = 1 - (1 - abx.r) ^ (1 / timestep)
  }
  
  # pi_Sr= probability of R transmitting to S (a proportion of pi_ssr if being colonised with S protects colonisation by R)
  pi_Sr = pi_ssr * (1 - bif)
  
  # For each day (first day should be filled)
  for(i in 2:nrow(patient.matrix)){
    
    # Get the previous row subset of the entire patient matrix which represents previous day or timestep
    prev_step = colo.matrix[i-1, ]
    # Get the indices which are already filled by a patient entering the ward for exclusion
    already_filled = which(!is.na(colo.matrix[i, ]))
    # Get the previous row subset of the entire abx matrix which represents previous day or timestep
    prev_abx = abx.matrix[i-1, ]
    # count if there are any R on the previous day
    r_num = length(which(prev_step == "R"))
    
    # probability for transmission of R to S
    infect_Sr = 1 - ((1 - pi_Sr) ^ (r_num / n.bed))
    # probability for transmission of R to ss
    infect_ssr = 1 - ((1 - pi_ssr) ^ (r_num / n.bed))
    
    if(abx.s + infect_Sr > 1){
      stop(paste("Error: abx.s + infect_Sr > 1 in scenario: those with broad antibiotics on previous day"))
    }
    
    if(repop.s + infect_ssr > 1){
      stop(paste("Error: repop.s + infect_ssr > 1 in scenario: those with no antibiotics on previous day"))
    }
    
    # First scenario: those with no antibiotics on previous day 
    ##########################################
    id_noabx = which(prev_abx == 0)
    
    ### Update R
    # Get the column indices which contain R in the previous day
    R = id_noabx[id_noabx %in% which(prev_step == "R")]
    # Remove column indices that already have a starting bacterial state filled in
    R = R[!(R %in% already_filled)]
    # if there is any R (number of R > 0) in the previous day
    if(length(R)){
      # Roll a random number for each R on the previous day
      roll = runif(length(R), 0, 1)
      # All column indices which roll < mu (decolonization parameter) are saved to fill in as S
      S_idx = R[roll < mu]
      # All the remaining column indices are saved to fill in as staying R the next day
      same_idx = R[roll >= mu]
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
      # Roll a random number for each S on the previous day for clearance
      roll= runif(length(S), 0, 1)
      
      r_idx = S[roll < infect_Sr]
      same_idx = S[roll >= infect_Sr]
      colo.matrix[i, r_idx] = "R"
      colo.matrix[i, same_idx] = "S"
    }
    
    ### Update ss
    ss = id_noabx[id_noabx %in% which(prev_step == "ss")]
    ss = ss[!(ss %in% already_filled)]
    if(length(ss)){
      roll = runif(length(ss), 0, 1)
      
      r_idx = ss[roll < infect_ssr]
      s_idx = ss[roll >= infect_ssr & roll < repop.s + infect_ssr]
      same_idx = ss[!(ss %in% c(r_idx, s_idx))]
      
      colo.matrix[i, s_idx] = "S"
      colo.matrix[i, r_idx] = "R"
      colo.matrix[i, same_idx] = "ss"
    }
    
    # Second scenario: those with narrow antibiotics on previous day 
    ##########################################
    id_narrowabx = which(prev_abx==1)
    
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
      # Roll a random number for each S on the previous day for clearance
      roll= runif(length(S), 0, 1)
      
      r_idx = S[roll < infect_Sr]
      ss_idx = S [roll >= infect_Sr & roll < abx.s + infect_Sr]
      same_idx= S[!(S %in% c(r_idx, ss_idx))]
      
      colo.matrix[i, ss_idx] = "ss"
      colo.matrix[i, r_idx] = "R"
      colo.matrix[i, same_idx] = "S"
    }
    
    ### Update ss
    ss = id_narrowabx[id_narrowabx %in%which(prev_step == "ss")]
    ss = ss[!(ss %in% already_filled)]
    if(length(ss)){
      roll = runif(length(ss), 0, 1)
      
      r_idx = ss[roll < infect_ssr]
      same_idx = ss[roll >= infect_ssr]
      
      colo.matrix[i, r_idx] = "R"
      colo.matrix[i, same_idx] = "ss"
    }
    
    # Third scenario: those with broad antibiotics on previous day 
    ##########################################
    id_broadabx = which(prev_abx == 2)
    
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
      ss_idx = R[roll < abx.r]
      same_idx = R[roll >= abx.r]
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
      # Roll a random number for each S on the previous day for clearance
      roll= runif(length(S), 0, 1)
      
      r_idx = S[roll < infect_Sr]
      ss_idx = S [roll >= infect_Sr & roll < abx.s + infect_Sr]
      same_idx= S[!(S %in% c(r_idx, ss_idx))]
      
      colo.matrix[i, ss_idx] = "ss"
      colo.matrix[i, r_idx] = "R"
      colo.matrix[i, same_idx] = "S"
    }
    
    ### Update ss
    ss = id_broadabx[id_broadabx %in% which(prev_step == "ss")]
    ss = ss[!(ss %in% already_filled)]
    if(length(ss)){
      roll = runif(length(ss), 0, 1)
      
      r_idx = ss[roll < infect_ssr]
      same_idx = ss[roll >= infect_ssr]
      
      colo.matrix[i, r_idx] = "R"
      colo.matrix[i, same_idx] = "ss"
    }
  }
  
  return(colo.matrix)
}

diff_prevalence <- function(n.bed, max.los, 
                            prop_R, prop_S, 
                            bif, pi_ssr, repop.s, mu, abx.s, abx.r,
                            p.infect, cum.r.1, p.r.day1, short_dur, long_dur){
  
  print(paste(c('parameters =', n.bed, max.los, 
                prop_R, prop_S, 
                bif, pi_ssr, repop.s, mu, abx.s, abx.r,
                p.infect, cum.r.1, p.r.day1, short_dur, long_dur)))
  
  old = Sys.time() # get start time
  
  timestep = 1
  n.day = 300
  iterations = 100
  
  iter_totalR = matrix(NA, nrow = n.day, ncol = iterations)
  for(iter in 1:iterations){
    matrixes = los.abx.table(n.bed = n.bed, n.day = n.day, max.los = max.los, 
                             p.infect = p.infect, p.r.day = p.r.day1, cum.r.1 = cum.r.1, 
                             meanDur = short_dur, timestep = timestep)
    patient.matrix = matrixes[[1]]
    if (all(!is.na(patient.matrix))) {print(paste('patient.matrix has no NA'))}
    abx.matrix = matrixes[[2]]
    if (all(!is.na(abx.matrix))) {print(paste('abx.matrix has no NA'))}
    los.array = summary.los(patient.matrix = patient.matrix)
    colo.matrix = colo.table(patient.matrix = patient.matrix, los = los.array, 
                             prop_R = prop_R, prop_S = prop_S)
    print(paste('colo.matrix is done'))
    colo_table_filled_iter = nextDay(patient.matrix = patient.matrix, los.array = los.array, 
                                     abx.matrix = abx.matrix, colo.matrix = colo.matrix, 
                                     bif = bif, pi_ssr = pi_ssr, repop.s = repop.s, mu = mu, 
                                     abx.s = abx.s, abx.r = abx.r, timestep = timestep)
    #Summary
    iter_totalR[, iter] = rowMeans(matrix(rowSums(colo_table_filled_iter == "R"), ncol = timestep, byrow=T))
    
    print(paste('colo.table is filled and results for short iter (', iter, ') =' , mean(iter_totalR[, iter][-c(1:150)]))) #mean of R for the last 150 days
  }
  # Discard first 150 runs as burn-in
  totalR_short = mean(rowSums(iter_totalR[151:nrow(iter_totalR), ,drop=FALSE]) / iterations / n.bed)
  
  
  iter_totalR = matrix(NA, nrow = n.day, ncol = iterations)
  for(iter in 1:iterations){
    matrixes = los.abx.table(n.bed = n.bed, n.day = n.day, max.los = max.los, 
                             p.infect = p.infect, p.r.day1 = p.r.day1, cum.r.1 = cum.r.1, 
                             meanDur = long_dur, timestep = timestep)
    patient.matrix = matrixes[[1]]
    abx.matrix = matrixes[[2]]
    los.array = summary.los(patient.matrix = patient.matrix)
    colo.matrix = colo.table(patient.matrix = patient.matrix, los = los.array, 
                             prop_R = prop_R, prop_S = prop_S)
    colo_table_filled_iter = nextDay(patient.matrix = patient.matrix, los.array = los.array, 
                                     abx.matrix = abx.matrix, colo.matrix = colo.matrix, 
                                     bif = bif, pi_ssr = pi_ssr, repop.s = repop.s, mu = mu, 
                                     abx.s = abx.s, abx.r = abx.r, timestep = timestep)
    
    #Summary
    iter_totalR[, iter] = rowMeans(matrix(rowSums(colo_table_filled_iter == "R"), ncol=timestep, byrow=T)) #mean of all the R for the time step/ sum of R for a day 
    
    print(paste('results for long iter (', iter, ') =' , mean(iter_totalR[, iter][-c(1:150)]))) #mean of R for the last 150 days
  }
  # Discard first 150 runs as burn-in
  totalR_long = mean(rowSums(iter_totalR[151:nrow(iter_totalR), ,drop=FALSE]) / iterations / n.bed)
  
  out = array(c(totalR_long, totalR_short, totalR_long - totalR_short))
  
  print(paste("overall mean difference", totalR_long - totalR_short))
  
  return(out)
  
}

prevalence <- function(n.bed, max.los, 
                       prop_R, prop_S, 
                       bif, pi_ssr, repop.s, mu, abx.s, abx.r,
                       p.infect, cum.r.1, p.r.day1, meanDur){
  
  old = Sys.time() # get start time
  
  timestep = 1
  iterations = 100
  n.day = 300
  
  iter_totalR = matrix(NA, nrow = n.day, ncol = iterations)
  
  for(iter in 1:iterations){
    matrixes = los.abx.table(n.bed = n.bed, n.day = n.day, max.los = max.los, 
                             p.infect = p.infect, p.r.day1 = p.r.day1, cum.r.1 = cum.r.1, 
                             meanDur= meanDur, timestep = timestep)
    patient.matrix = matrixes[[1]]
    abx.matrix = matrixes[[2]]
    los.array = summary.los(patient.matrix = patient.matrix)
    colo.matrix = colo.table(patient.matrix = patient.matrix, los = los.array, 
                             prop_R = prop_R,prop_S = prop_S)
    
    colo_table_filled_iter = nextDay(patient.matrix = patient.matrix, los.array = los.array, 
                                     abx.matrix = abx.matrix, colo.matrix = colo.matrix, 
                                     bif = bif, pi_ssr = pi_ssr, repop.s = repop.s, mu = mu, 
                                     abx.s = abx.s, abx.r = abx.r,timestep = timestep)

    #Summary
    iter_totalR[, iter] = rowMeans(matrix(rowSums(colo_table_filled_iter == "R"), ncol = timestep, byrow = T)) #rowMeans(number of days by timestep) gives mean number of R per day
  }
  # Discard first 150 runs as burn-in
  totalR = mean(rowSums(iter_totalR[151:nrow(iter_totalR), , drop = FALSE]) / iterations / n.bed) #mean(avg R per bed per day)
  
  #print(paste("totalR_long", totalR_long, "totalR_short", totalR_short))
  # print elapsed time
  new = Sys.time() - old # calculate difference
  print(new) # print in nice format
  return(totalR)
}

parameters_prevalence_simple<- c("n.bed", "max.los", 
                                 "prop_R", "prop_S", 
                                 "bif", "pi_ssr", "repop.s", "mu", 
                                 "abx.s", "abx.r", "p.infect", "cum.r.1", 'p.r.day1', "meanDur")

parameters_diff_prevalence_simple<- c("n.bed", "max.los", 
                                      "prop_R", "prop_S", 
                                      "bif", "pi_ssr", "repop.s", "mu", 
                                      "abx.s", "abx.r", "p.infect", "cum.r.1", 'p.r.day1', "short_dur", "long_dur")


