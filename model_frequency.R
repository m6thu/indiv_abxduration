source('msm_util_rtnorm.R')
source('los_abx_matrix.R')

r_beta = rbeta(5000, 0.2, 2) #distribution of R in the carriers
r_beta_norm = (r_beta-min(r_beta)) / (max(r_beta) - min(r_beta)) #normalised from 0 to 1

# shape 1 and shape 2  based on rene gut data 
# d = read.csv('gutdata/Ini_CTXm_copies_qPCR.csv')
# hist(d$ini_CTXm_copies)
# hist(rbeta(10000, 0.1, 2))
# hist(r_beta_norm*exp(16))

colo.table <- function(patient.matrix, los.array, total_prop, prop_R, r_thres, K){
  
  n.day = nrow(patient.matrix)
  n.bed = ncol(patient.matrix)
  
  number_of_patients = dim(los.array)[2]
  
  #capacity for enterobacteriaceae growth (log)
  #log of the capacity is normal in distribution from rene 
  #d=read.csv('gutdata/Ini_CTXm_copies_qPCR.csv')
  #hist(d$ini_16S_log)
  total_capacity = rnorm(number_of_patients, mean = K) #in log 
  total_capacity_matrix = matrix(rep(total_capacity, los.array[2,]), byrow = F, ncol = ncol(patient.matrix)) #in log 
  
  #existing population 
  #existing population mean is proportion (total_prop) of total capacity instead of a proportion of the 
  #distribution of the total capacity so that some starts at full capacity while others are not 
  #which is similar to model 2 where there is sr 
  total_existing_mean = log(total_prop * exp(total_capacity)) #in log 
  total_existing = rnorm(number_of_patients, mean = total_existing_mean) #total number of Enterobacteriaceae is a proportion of the capacity (log)
  total_existing [which(total_existing >= total_capacity)] = total_capacity [which(total_existing >= total_capacity)] #those exceeding total capacity will be given their own full capacity
  
  #amount of S and R carried 
  morethan_rthres = which(total_existing > r_thres)  #pick patients who have existing Enterobacteriaceae more than r_thres (log) 
  r.id = sample(morethan_rthres, size= floor(prop_R * number_of_patients)) #potential bug here- morethan_rthres may be less than size of r_id
  if (length(r.id) > 0) {s.id = (1:number_of_patients) [-r.id]} else {s.id = (1:number_of_patients)}
  
  r_bact = rep(NA, number_of_patients)
  # R individuals carry R 
  for (ind in r.id){  #in log 
    #expand r_beta to each individual's existing population
    r_list_full = r_beta_norm * total_existing [ind] 
    r_bact[ind] = sample(r_list_full[which(r_list_full >= r_thres)], 1)#list of r amounts to sample from for R
  }
  
  # S individuals also carry some R (similar to Sr, sr, sR in model 2)
  for (ind in s.id){  #in log 
    #expand r_beta to each individual's existing population
    r_list_full = r_beta_norm * total_existing [ind]
    r_bact[ind] = sample(r_list_full[which(r_list_full < r_thres)], 1)#list of r amounts to sample from for S
  }
  
  r_bact[which(r_bact<0)] = -Inf #those with less than 1 bacteria round to 0 (log)
  
  #check below equal to prop_R
  #sum(r_bact>r_thres)/number_of_patients
  #check shape 
  #hist(exp(r_bact))
  s_bact = exp(total_existing) - exp(r_bact) #abs
  s_bat_abs = ifelse(s_bact<0, 0, s_bact)
  s_bact = log(s_bat_abs) #total number of S for each patient (log)
  
  S_Bactlevelstart = matrix(NA, n.day, n.bed)
  R_Bactlevelstart = matrix(NA, n.day, n.bed)
  
  # pad with NAs
  end_idx = 1
  for(i in 1:number_of_patients){
    S_Bactlevelstart[end_idx:(end_idx + los.array[2, i] - 1)] = c(s_bact[i], rep(NA, los.array[2, i]-1))
    R_Bactlevelstart[end_idx:(end_idx + los.array[2, i] - 1)] = c(r_bact[i], rep(NA, los.array[2, i]-1))
    end_idx = end_idx + los.array[2, i]
  }
  
  return(list(S_Bactlevelstart, R_Bactlevelstart, total_capacity_matrix)) # in log
}

# Update values for every day (define function)
nextDay <- function(patient.matrix, los.array, abx.matrix, colo.matrix, 
                    pi_ssr, total_prop,  K, r_growth, r_thres, r_trans, s_growth,
                    abx.s, abx.r, timestep){
  
  if (timestep > 1) {
    pi_ssr = 1 - (1 - pi_ssr) ^ (1 / timestep)
  }
  
  S_table = colo.matrix[[1]] #in log
  R_table = colo.matrix[[2]] #in log
  
  #capacity matrix for enterobacteriaceae growth (log)
  total_capacity = colo.matrix[[3]]
  
  # For each day (first day should be filled)
  for(i in 2:nrow(patient.matrix)){
    
    # calculate how many people has R above r_thres (log)
    r_num = sum(R_table[i-1,] >= r_thres)
    n.bed = ncol(patient.matrix)
    
    if (is.na(r_num)) {r_num = 0}
    
    # from number of r, calculate probability of transmission
    prop_r = 1 - ((1 - pi_ssr) ^ (r_num / n.bed)) 
    
    ###### Convert all log scale parameters into normal scale for addition, then convert back to log
    #for each person:
    for(j in 1:ncol(patient.matrix)){
      
      #print(c(i,j))
      
      if(is.na(R_table[i, j])){ # pick any; S and R should be filled in same slots
        # calculate effect of R logistic bacteria growth (abs) - only R growth if on antibiotics
        select_pa = as.numeric(abx.matrix[i-1, j] != 0)
        R_grow = select_pa * r_growth * exp(R_table[i-1, j]) * (1 - ((exp(R_table[i-1, j]) + exp(S_table[i-1, j]))/exp(total_capacity[i, j]))) 
        # abx killing if abx.matrix is r abx (== 2) (abs)
        R_abx = -(abx.matrix[i-1, j] == 2) * abx.r * exp(R_table[i-1, j])
        # add effects to current table (abs first because log of a negative number is NaN)
        R_table[i, j] = exp(R_table[i-1, j]) + R_grow + R_abx
      } else {
        R_table[i, j] = exp(R_table[i, j]) #if filled, change to abs
      }
      
      if(is.na(S_table[i, j])){ 
        # calculate effect of S logistic bacteria growth (in absolute numbers)
        S_grow = s_growth * exp(S_table[i-1, j]) * (1 - ((exp(R_table[i-1, j]) + exp(S_table[i-1, j]))/exp(total_capacity[i, j])))
        # calculate killing effect of antibiotics R and antibiotics S (abs)
        S_abx = -(abx.matrix[i-1, j] >= 1) * abx.s * exp(S_table[i-1, j])
        # apply effects (abs first because log of a negative number is NaN)
        S_table[i, j] = exp(S_table[i-1, j]) + S_grow + S_abx
      } else { 
        S_table[i, j] = exp(S_table[i, j]) #if filled, change to abs
      }
      
      #transmission of R if roll pass prob check 
      roll = runif(1, 0, 1) # roll for transmission
      R_trans = exp(r_trans) * (roll < prop_r) # abs 
      
      # transmission and trim range 
      ### transmission only happens if R and S have not exceeded total capacity and S is not fully occupying the capacity
      if(S_table[i, j] + R_table[i, j] > exp(total_capacity[i, j])){ ## if existing S and R already exceed total capacity (abs)
        
        #if (R_table[i, j]<1) { ###  if the whole capacity is occupied by S, transmission happens with R_trans
        #   R_table[i, j] = R_table[i, j] + R_trans 
        #   # this will make S+R to exceed total capacity but next day S and R will be reduced proportionally
        
        #} else { ### no transmission if S and R more than total_capacity
        S_table[i, j] = S_table[i, j] / (S_table[i, j] + R_table[i, j]) * exp(total_capacity[i, j]) #S and R each given their proportion of the total capacity (abs)
        R_table[i, j] = R_table[i, j] / (S_table[i, j] + R_table[i, j]) * exp(total_capacity[i, j]) 
        #natural attrition of S and R take place according to their density if capacity exceeded
        
        #}
        
      } else { #transmission if S and R less than total_capacity
        
        if (S_table[i, j] + R_table[i, j] + R_trans > exp(total_capacity[i, j])) {
          ### transmission if S and R less than total_capacity but with transmission exceeds total_capacity
          S_table[i, j] = S_table[i, j] / (S_table[i, j] + R_table[i, j] + R_trans) * exp(total_capacity[i, j]) #S and R+R_trans each given their proportion of the total capacity (abs)
          R_table[i, j] = (R_table[i, j] + R_trans) / (S_table[i, j] + R_table[i, j] + R_trans) * exp(total_capacity[i, j]) 
          
        } else {
          ### transmission if sum of S and R and transmitted R are less than the total capacity
          #transmission happens with full amount of R_trans
          R_table[i, j] = R_table[i, j] + R_trans 
        }
      }
      
      # covert to log 
      if (S_table[i, j]<0) {S_table[i, j] = log(0)} else {S_table[i, j] =log(S_table[i, j])}
      if (R_table[i, j]<0) {R_table[i, j] = log(0)} else {R_table[i, j] =log(R_table[i, j])}
    }
  }
  return(list(S_table, R_table))
}

diff_prevalence <- function(n.bed, max.los, p.infect, cum.r.1, p.r.day1,
                            K, total_prop,  prop_R, pi_ssr, 
                            r_trans, r_growth, r_thres, s_growth,
                            abx.s, abx.r, short_dur,long_dur){
  
  old = Sys.time() # get start time
  print(old)
  # DEBUG
  print('parameter values')
  print(paste(n.bed, max.los, p.infect, cum.r.1, p.r.day1,
              K, total_prop,  prop_R, pi_ssr, r_trans, r_growth, r_thres, s_growth,
              abx.s, abx.r, short_dur,long_dur))
  
  timestep = 1
  n.day = 300
  iterations=100
  
  #iter_totalR.no = matrix(NA, nrow = n.day, ncol = iterations)
  iter_totalR.thres = matrix(NA, nrow = n.day, ncol = iterations)
  
  for(iter in 1:iterations){
    matrixes = los.abx.table(n.bed=n.bed, n.day=n.day, max.los=max.los, 
                             p.infect=p.infect, p.r.day1=p.r.day1, cum.r.1=cum.r.1, 
                             meanDur= short_dur, timestep=timestep)
    patient.matrix=matrixes[[1]]
    abx.matrix=matrixes[[2]]
    los.array = summary.los(patient.matrix=patient.matrix)
    colo.matrix = colo.table(patient.matrix=patient.matrix, los.array=los.array, total_prop=total_prop, prop_R=prop_R, r_thres=r_thres, K=K)
    
    colo.matrix_filled_iter = nextDay(patient.matrix=patient.matrix, los.array=los.array, abx.matrix=abx.matrix, colo.matrix=colo.matrix, 
                                      pi_ssr=pi_ssr, total_prop=total_prop, K=K, 
                                      r_growth=r_growth, r_thres=r_thres,r_trans=r_trans, s_growth=s_growth,
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
  totalR_thres_short = mean(rowSums(iter_totalR.thres[151:nrow(iter_totalR.thres),, drop=FALSE])/iterations/n.bed)
  
  #iter_totalR.no = matrix(NA, nrow = n.day, ncol = iterations)
  iter_totalR.thres = matrix(NA, nrow = n.day, ncol = iterations)
  
  for(iter in 1:iterations){
    matrixes = los.abx.table(n.bed=n.bed, n.day=n.day, max.los=max.los, 
                             p.infect=p.infect, p.r.day1=p.r.day1, cum.r.1=cum.r.1, 
                             meanDur= long_dur, timestep=timestep)
    patient.matrix=matrixes[[1]]
    abx.matrix=matrixes[[2]]
    los.array = summary.los(patient.matrix=patient.matrix)
    colo.matrix = colo.table(patient.matrix=patient.matrix, los.array=los.array, total_prop=total_prop, prop_R=prop_R, r_thres=r_thres, K=K)
    
    colo.matrix_filled_iter = nextDay(patient.matrix=patient.matrix, los.array=los.array, abx.matrix=abx.matrix, colo.matrix=colo.matrix, 
                                      pi_ssr=pi_ssr, total_prop=total_prop, K=K, r_trans=r_trans,
                                      r_growth=r_growth, r_thres=r_thres, s_growth=s_growth,
                                      abx.s=abx.s, abx.r=abx.r, timestep=timestep)
    
    #Summary 
    #for total units of R bacteria on a day
    df.R = data.frame(colo.matrix_filled_iter[[2]])
    #iter_totalR.no[, iter] = rowMeans(matrix(rowSums(df.R), ncol=timestep, byrow=T))
    #for number of people who reached R threshold on a day
    iter_totalR.thres[, iter] = rowMeans(matrix(rowSums(df.R >= r_thres), ncol=timestep, byrow=T))
    #print("end iteration loop")
  }
  #totalR_no_long = mean(rowSums(iter_totalR.no[151:nrow(iter_totalR.no),, drop=FALSE])/iterations/n.bed)
  totalR_thres_long = mean(rowSums(iter_totalR.thres[151:nrow(iter_totalR.thres),, drop=FALSE])/iterations/n.bed)
  
  print('calculation complete')
  # print elapsed time
  new = Sys.time() - old # calculate difference
  print(new) # print in nice format
  
  print(c(totalR_thres_long, totalR_thres_short, totalR_thres_long-totalR_thres_short))
  
  return(array(c(totalR_thres_long, totalR_thres_short, totalR_thres_long-totalR_thres_short)))
}

prevalence <- function(n.bed, max.los, p.infect, cum.r.1, p.r.day1,
                       K, total_prop,  prop_R, pi_ssr, r_trans, r_growth, r_thres, s_growth,
                       abx.s, abx.r, meanDur){
  
  old = Sys.time() # get start time
  # DEBUG
  print(paste(n.bed, max.los, p.infect, cum.r.1, p.r.day1,
              K, total_prop,  prop_R, pi_ssr, r_trans, r_growth, r_thres, s_growth,
              abx.s, abx.r, meanDur))
  
  timestep = 1
  n.day = 300
  iterations= 100
  
  iter_totalR.thres = matrix(NA, nrow = n.day, ncol = iterations)
  
  for(iter in 1:iterations){
    
    matrixes = los.abx.table(n.bed=n.bed, n.day=n.day, max.los=max.los, 
                             p.infect=p.infect, p.r.day1=p.r.day1, cum.r.1=cum.r.1, 
                             meanDur= meanDur, timestep=timestep)
    patient.matrix=matrixes[[1]]
    abx.matrix=matrixes[[2]]
    los.array = summary.los(patient.matrix=patient.matrix)
    colo.matrix = colo.table(patient.matrix=patient.matrix, los.array=los.array, total_prop=total_prop, prop_R=prop_R, r_thres=r_thres, K=K)
    
    colo.matrix_filled_iter = nextDay(patient.matrix=patient.matrix, los.array=los.array, abx.matrix=abx.matrix, colo.matrix=colo.matrix, 
                                      pi_ssr=pi_ssr, total_prop=total_prop, K=K, r_growth=r_growth,r_thres=r_thres, r_trans=r_trans,s_growth=s_growth,
                                      abx.s=abx.s, abx.r=abx.r, timestep=timestep)
    # Summary
    df.R = data.frame(colo.matrix_filled_iter[[2]])
    #for number of people who reached R threshold on a day
    iter_totalR.thres[, iter]=rowMeans(matrix(rowSums(df.R >= r_thres), ncol=timestep, byrow=T))
    #print("end iteration loop")
  }
  
  totalR_thres = mean(rowSums(iter_totalR.thres[151:nrow(iter_totalR.thres), , drop=FALSE])/iterations/n.bed)
  
  # print elapsed time
  new = Sys.time() - old # calculate difference
  print(new) # print in nice format
  
  return(totalR_thres)
}

parameters_prevalence_freq <- c("n.bed", "max.los", "p.infect", "cum.r.1", "p.r.day1", 
                                "K", "total_prop", "prop_R",
                                "pi_ssr", "r_trans", "r_growth", 'r_thres', 's_growth',
                                "abx.s", "abx.r",
                                "meanDur")

parameters_diff_prevalence_freq <- c("n.bed", "max.los", "p.infect", "cum.r.1", "p.r.day1", 
                                     "K", "total_prop", "prop_R",
                                     "pi_ssr", "r_trans", "r_growth", 'r_thres','s_growth',
                                     "abx.s", "abx.r",
                                     "short_dur", "long_dur")

