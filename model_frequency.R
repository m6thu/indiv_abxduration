source('msm_util_rtnorm.R')
source('los_abx_matrix.R')

colo.table <- function(patient.matrix, los.array, total_prop, prop_R, r_mean, K){
  
  n.day = nrow(patient.matrix)
  n.bed = ncol(patient.matrix)
  
  number_of_patients = dim(los.array)[2]
  
  #total capacity for enterobacteriaceae growth (log)
  total_capacity = rnorm(number_of_patients, K)
  total_capacity_matrix= matrix(rep(total_capacity, los.array[2,]), byrow = F, ncol = ncol(patient.matrix))
  
  #existing population 
  total_prop= rnorm(number_of_patients,total_prop)
  total_prop.norm= (total_prop-min(total_prop))/(max(total_prop)-min(total_prop)) #normalised proportion 
  total_existing= log(total_prop.norm*exp(total_capacity)) #total number of Enterobacteriaceae is a proportion of the total capacity (log)
  
  prop_R_id = sample(1:number_of_patients, prop_R*number_of_patients) #patients who carry resistant organisms
  prop.r.ent= rtnorm(number_of_patients, mean=r_mean, lower=0)
  prop.r.ent.norm= (prop.r.ent-min(prop.r.ent))/(max(prop.r.ent)-min(prop.r.ent))#proportion of R in total population of Enterobacteriaceae (normalised)
  r_transmitted= log(prop.r.ent.norm*exp(total_existing))
  r_transmit_matrix= matrix(rep(r_transmitted, los.array[2,]), byrow = F, ncol = ncol(patient.matrix))
  
  r_bact = rep(log(0), number_of_patients)
  r_bact[prop_R_id]= log(prop.r.ent.norm[prop_R_id]*exp(total_existing[prop_R_id])) #total number of resistant Enterobacteriaceae for each patient (log)
  s_bact = log(exp(total_existing) - exp(r_bact)) #total number of sensitive Enterobacteriaceae for each patient (log)
  
  S_Bactlevelstart = matrix(NA, n.day, n.bed)
  R_Bactlevelstart = matrix(NA, n.day, n.bed)
  
  # pad with NAs
  end_idx = 1
  for(i in 1:number_of_patients){
    S_Bactlevelstart[end_idx:(end_idx + los.array[2, i] - 1)] = c(s_bact[i], rep(NA, los.array[2, i]-1))
    R_Bactlevelstart[end_idx:(end_idx + los.array[2, i] - 1)] = c(r_bact[i], rep(NA, los.array[2, i]-1))
    end_idx = end_idx + los.array[2, i]
  }
  
  return(list(S_Bactlevelstart, R_Bactlevelstart, total_capacity_matrix, r_transmit_matrix)) # in log
}

# Update values for every day (define function)
nextDay <- function(patient.matrix, los.array, abx.matrix, colo.matrix, 
                    pi_ssr, total_prop,  K, r_mean, r_growth, r_thres, s_growth,
                    abx.s, abx.r, timestep){
  
  if (abx.r<0.05) { #if abx.r is ineffective in scenario B - resistance = CRE 
    abx.r.aginst.s =abx.s #remains effective for s  
    abx.r.aginst.r =abx.r #ineffective for r 
  } else { # if abx.r is effective in scenario B - resistance = ESBL
    abx.r.aginst.s =abx.r #effective for both s and r 
    abx.r.aginst.r =abx.r
  }
  
  # K: loading capacity
  pi_ssr = 1-(1-pi_ssr)^(1/timestep)
  
  S_table = colo.matrix[[1]] #in log
  R_table = colo.matrix[[2]] #in log
  
  #total capacity matrix for enterobacteriaceae growth (log)
  total_capacity=colo.matrix[[3]]
  
  #amount of R transferred 
  r_trans=colo.matrix[[4]]
  
  #threshold of Enterobacteriaceae before the patient can transmit 
  r_thres_matrix= log(r_thres* exp(total_capacity))
  
  # For each day (first day should be filled)
  for(i in 2:nrow(patient.matrix)){
    # calculate how many people has R above 0 (log)
    r_num = sum(R_table[i-1,] >= r_thres_matrix[i-1,])
    # from number of r, calculate probability of transmission
    prop_r = 1-((1-pi_ssr)^r_num) 
    
    ###### Convert all log scale parameters into normal scale for addition, then convert back to log
    #for each person:
    for(j in 1:ncol(patient.matrix)){
      if(is.na(R_table[i, j])){ # pick any; S and R should be filled in same slots
        
        # roll for transmission
        roll = runif(1, 0, 1)
        # calculate effect of R logistic bacteria growth (abs)
        R_grow = r_growth*exp(R_table[i-1, j])*(1 - ((exp(R_table[i-1, j]) + exp(S_table[i-1, j]))/exp(total_capacity[i, j]))) 
        #add effect of transmission if roll pass prob check and if previous R level is 0 (abs)
        R_trans = exp(r_trans[i,j])*(roll < prop_r)# & !R_table[i-1, j])
        # add effect of abx death if abx.matrix is r abx (== 2) (abs)
        R_abx = -(abx.matrix[i-1, j] == 2)*abx.r.aginst.r*R_grow
        # apply effects to current table (abs first because log of a negative number is NaN)
        R_table[i, j] = exp(R_table[i-1, j]) + R_grow + R_trans + R_abx
        
        # if (R_grow>0){
        #   print(paste0('capacity',exp(total_capacity[i, j])))
        #   print(paste0('old R',exp(R_table[i-1, j])))
        #   print(paste0('R_grow',R_grow))
        #   print(paste0('old S',exp(S_table[i-1, j])))
        #   print(paste0('new R',exp(R_table[i, j])))
        # }
        
        # trim if value goes beyond range
        if(R_table[i, j] > exp(total_capacity[i, j])){
          R_table[i, j] = exp(total_capacity[i, j]) #abs
        }
        if(R_table[i, j] < 1){ #if less than 1 bacteria
          R_table[i, j] = 0 #abs
        }
        R_table[i, j] = log(R_table[i, j]) #log
      }
      
      if(is.na(S_table[i, j])){ 
        # calculate effect of S logistic bacteria growth (in absolute numbers)
        S_grow = s_growth*exp(S_table[i-1, j])*(1 - ((exp(R_table[i-1, j]) + exp(S_table[i-1, j]))/exp(total_capacity[i, j])))
        # calculate effect of death antibiotics R and effect of death by abx S (abs)
        S_abx_s = -(abx.matrix[i-1, j] == 1)*abx.s*S_grow
        S_abx_r = -(abx.matrix[i-1, j] == 2)*abx.r.aginst.s*S_grow
        # apply effects (abs first because log of a negative number is NaN)
        S_table[i, j] = exp(S_table[i-1, j]) + S_grow + S_abx_s + S_abx_r
        # trim range
        if(S_table[i, j] > exp(total_capacity[i, j])){
          S_table[i, j] = exp(total_capacity[i, j])
        }
        if(S_table[i, j] < 1){ #if less than 1 bacteria
          S_table[i, j] = 0
        }
        S_table[i, j] = log(S_table[i, j]) #log
      }
    }
  }
  
  return(list(S_table, R_table))
}

diff_prevalence <- function(n.bed, max.los, p.infect, cum.r.1, p.r.day1,
                            K, total_prop,  prop_R, pi_ssr, 
                            r_mean, r_growth, r_thres, s_growth,
                            abx.s, abx.r, short_dur,long_dur){
  
  old = Sys.time() # get start time
  # DEBUG
  print(paste(n.bed, max.los, p.infect, cum.r.1, p.r.day1,
              K, total_prop,  prop_R, pi_ssr, r_mean, r_growth, r_thres, s_growth,
              abx.s, abx.r, short_dur,long_dur))
  
  timestep = 1
  n.day = 300
  
  if (abx.r> 1){ #scenario A
    iterations= 100
  } else { #scenario B
    iterations= 50
  }
 
  
  #iter_totalR.no = matrix(NA, nrow = n.day, ncol = iterations)
  iter_totalR.thres = matrix(NA, nrow = n.day, ncol = iterations)
  
  for(iter in 1:iterations){
    
    matrixes = los.abx.table(n.bed=n.bed, n.day=n.day, max.los=max.los, 
                             p.infect=p.infect, p.r.day1=p.r.day1, cum.r.1=cum.r.1, 
                             meanDur= short_dur, timestep=timestep)
    patient.matrix=matrixes[[1]]
    abx.matrix=matrixes[[2]]
    los.array = summary.los(patient.matrix=patient.matrix)
    colo.matrix = colo.table(patient.matrix=patient.matrix, los.array=los.array, total_prop=total_prop, prop_R=prop_R, r_mean=r_mean, K=K)
    
    colo.matrix_filled_iter = nextDay(patient.matrix=patient.matrix, los.array=los.array, abx.matrix=abx.matrix, colo.matrix=colo.matrix, 
                                      pi_ssr=pi_ssr, total_prop=total_prop, K=K, 
                                      r_mean=r_mean, r_growth=r_growth, r_thres=r_thres, s_growth=s_growth,
                                      abx.s=abx.s, abx.r=abx.r, timestep=timestep)
    # Summary - timestep by bed in absolute numbers
    df.R = data.frame(colo.matrix_filled_iter[[2]])
    r_thres_matrix= data.frame(colo.matrix[[4]])
    # for number of people who reached R threshold on a day
    ##   sum of number of people per timestep that reach threshold 
    ##   make a matrix of sum of people per day (days by timestep)
    ##   daily means of number of people who reached R threshold 
    iter_totalR.thres[, iter]=rowMeans(matrix(rowSums(df.R >= r_thres_matrix), ncol=timestep, byrow=T))
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
    colo.matrix = colo.table(patient.matrix=patient.matrix, los.array=los.array, total_prop=total_prop, prop_R=prop_R, r_mean=r_mean, K=K)
    
    colo.matrix_filled_iter = nextDay(patient.matrix=patient.matrix, los.array=los.array, abx.matrix=abx.matrix, colo.matrix=colo.matrix, 
                                      pi_ssr=pi_ssr, total_prop=total_prop, K=K, 
                                      r_mean=r_mean, r_growth=r_growth, r_thres=r_thres, s_growth=s_growth,
                                      abx.s=abx.s, abx.r=abx.r, timestep=timestep)
    
    #Summary 
    #for total units of R bacteria on a day
    df.R = data.frame(colo.matrix_filled_iter[[2]])
    r_thres_matrix= data.frame(colo.matrix[[4]])
    #iter_totalR.no[, iter] = rowMeans(matrix(rowSums(df.R), ncol=timestep, byrow=T))
    #for number of people who reached R threshold on a day
    iter_totalR.thres[, iter] = rowMeans(matrix(rowSums(df.R >= r_thres_matrix), ncol=timestep, byrow=T))
    #print("end iteration loop")
  }
  #totalR_no_long = mean(rowSums(iter_totalR.no[151:nrow(iter_totalR.no),, drop=FALSE])/iterations/n.bed)
  totalR_thres_long = mean(rowSums(iter_totalR.thres[151:nrow(iter_totalR.thres),, drop=FALSE])/iterations/n.bed)
  
  # print elapsed time
  new = Sys.time() - old # calculate difference
  print(new) # print in nice format
  
  return(array(c(totalR_thres_long, totalR_thres_short, totalR_thres_long-totalR_thres_short)))
}

prevalence <- function(n.bed, max.los, p.infect, cum.r.1, p.r.day1,
                       K, total_prop,  prop_R, pi_ssr, r_mean, r_growth, r_thres, s_growth,
                       abx.s, abx.r, meanDur){
  
  old = Sys.time() # get start time
  # DEBUG
  print(paste(n.bed, max.los, p.infect, cum.r.1, p.r.day1,
              K, total_prop,  prop_R, pi_ssr, r_mean, r_growth, r_thres, s_growth,
              abx.s, abx.r, meanDur))
  
  timestep = 1
  n.day = 300
  
  if (abx.r > 1){ #scenario A
    iterations= 100
  } else { #scenario B
    iterations= 50
  }
  
  iter_totalR.thres = matrix(NA, nrow = n.day, ncol = iterations)
  
  for(iter in 1:iterations){
    
    matrixes = los.abx.table(n.bed=n.bed, n.day=n.day, max.los=max.los, 
                             p.infect=p.infect, p.r.day1=p.r.day1, cum.r.1=cum.r.1, 
                             meanDur= meanDur, timestep=timestep)
    patient.matrix=matrixes[[1]]
    abx.matrix=matrixes[[2]]
    los.array = summary.los(patient.matrix=patient.matrix)
    colo.matrix = colo.table(patient.matrix=patient.matrix, los.array=los.array, total_prop=total_prop, prop_R=prop_R,r_mean = r_mean,K=K)
    
    colo.matrix_filled_iter = nextDay(patient.matrix=patient.matrix, los.array=los.array, abx.matrix=abx.matrix, colo.matrix=colo.matrix, 
                                      pi_ssr=pi_ssr, total_prop=total_prop, K=K, r_mean=r_mean, r_growth=r_growth,r_thres=r_thres, s_growth=s_growth,
                                      abx.s=abx.s, abx.r=abx.r, timestep=timestep)
    # Summary
    df.R = data.frame(colo.matrix_filled_iter[[2]])
    r_thres_matrix= data.frame(colo.matrix[[4]])
    #for number of people who reached R threshold on a day
    iter_totalR.thres[, iter]=rowMeans(matrix(rowSums(df.R >= r_thres_matrix), ncol=timestep, byrow=T))
    #print("end iteration loop")
  }
  
  totalR_thres = mean(rowSums(iter_totalR.thres[151:nrow(iter_totalR.thres),, drop=FALSE])/iterations/n.bed)
  
  # print elapsed time
  new = Sys.time() - old # calculate difference
  print(new) # print in nice format
  
  return(totalR_thres)
}

parameters_prevalence_freq <- c("n.bed", "max.los", "p.infect", "cum.r.1", "p.r.day1", 
                                "K", "total_prop", "prop_R",
                                "pi_ssr", "r_mean", "r_growth", 'r_thres', 's_growth',
                                "abx.s", "abx.r",
                                "meanDur")

parameters_diff_prevalence_freq <- c("n.bed", "max.los", "p.infect", "cum.r.1", "p.r.day1", 
                                     "K", "total_prop", "prop_R",
                                     "pi_ssr", "r_mean", "r_growth", 'r_thres','s_growth',
                                     "abx.s", "abx.r",
                                     "short_dur", "long_dur")

