#(1 for short duration and 1 for long duration) 
#simulate inpatients with various lengths of stay
#allocate various duration of antibiotics for each patient

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

# abx.table <- function(patient.matrix, los.array, p.s, p.r.day1, p.r.dayafter, 
#                       meanDur.s, meanDur.r, sdDur, timestep=1){
#     
#     stopifnot(p.s+p.r.day1 <= 1)
#     
#     #generate antibiotic use table
#     #number of days of s antibiotic is randomly drawn from a truncated normal distribution
#     abx_days.s <- round(rtnorm(ncol(los.array), mean=meanDur.s*timestep, sd=sdDur*timestep, lower=0))
#     #number of days of r antibiotic is drawn from the distribution of accumulated probability
#     abx_days.r <- p.r.dayafter
#     # Unit test - check distribution of abx distribution
#     # hist(abx_days, breaks=20)
#     # Unit test - compare cases that will enter padding if-else
#     # abx_days > los.array[2, ]
#     
#     # abx.matrix should be same size as patient.matrix
#     abx.matrix <- matrix(NA, nrow=nrow(patient.matrix), ncol=ncol(patient.matrix))
#     idx_end <- 1
#     # For each patient
#     for(i in 1:ncol(los.array)){
#         
#         # maxiumum number of days for a particular patient from los vector
#         max_days <- los.array[2, i]
#         # number of abx days for that particular patient from generated number
#         abx_s <- abx_days.s[i]
#         abx_r <- abx_days.r[i]
#         # Initial treatment value derived from probability
#         rand <- runif(1, 0, 1)
#         if (rand < p.s){
#             if(abx_person > max_days){
#                 # if the number of generated abx is longer than the los of that person
#                 # have them take abx everyday for their stay
#                 abx.matrix[idx_end:(idx_end+max_days-1)] <- rep(1, max_days)
#             }else{
#                 # else take abx for abx days and pad to fit max_days
#                 abx.matrix[idx_end:(idx_end+max_days-1)] <- c(rep(1, abx_person), rep(0, max_days-abx_person))
#             }
#         }else{
#             # no abx taken for that person
#             abx.matrix[idx_end:(idx_end+max_days-1)] <- rep(0, max_days)
#         }
#         # move starting position to end of previous patient
#         idx_end <- idx_end+max_days
#     }
#     
#     return(abx.matrix)
# }

abx.table <- function(patient.matrix, los.array, p.s, p.r.day1, p.r.dayafter,
                      meanDur.s, meanDur.r, sdDur, timestep=1){

    n.day <- nrow(patient.matrix)
    n.bed <- ncol(patient.matrix)
    
    #generate antibiotic use table
    ## antibiotics to treat sensitive organisms 
    matrix_DuraDay <- matrix(NA, nrow=n.day, ncol=n.bed)
    #matrix_DuraDay is the matrix with number of days of antibiotics.s for every patient in patient.matrix
    
    get1stdayofstay <- function(patientnumber,bednumber, bedoccmat){
        match(patientnumber, bedoccmat[,bednumber])
    }
    #returns first day of stay for patientnumber in bednumber, or NA if patient not there
    
    #     getlos<-function(patientnumber,bednumber, bedoccmat){
    #         sum(bedoccmat[,bednumber]==patientnumber)
    #     }
    #returns los of patient patientnumber in bed bedumber, in bed occupancy matrix, 
    #0 if patient not found
    
    for (i in 1:max(patient.matrix)){
        for (j in 1:n.bed){
            matrix_DuraDay[get1stdayofstay(i,j,patient.matrix), j] <- abs(round(rnorm(1, mean=meanDur.s*timestep, sd=sdDur*timestep)))
        }
    }
    #number of days of antibiotic.s is randomly drawn from a normal dist
    
    for (i in 2:n.day){
        for (j in 1:n.bed){
            if(is.na(matrix_DuraDay[i,j])){
                matrix_DuraDay[i, j] <- matrix_DuraDay[i-1,j]
            }else{
                matrix_DuraDay[i, j] <- matrix_DuraDay[i, j]
            }
        }
    }
    #Fill the matrix for antibiotics.s with same dimension as patient.matrix
    
    matrix_AtbTrt <- matrix(NA, nrow=n.day, ncol=n.bed) 
    # #matrix_AtbTrt to count the cummulative length of stay for treated patients in the same dimension as patient.matrix
    
    for (i in 1:max(patient.matrix)){
        for (j in 1:n.bed){
            rand <- runif(1,0,1)
            if (rand < p.s) {
                matrix_AtbTrt[get1stdayofstay(i,j,patient.matrix), j] <-  1
            } else {
                matrix_AtbTrt[get1stdayofstay(i,j,patient.matrix), j] <-  0
            }
        }
    }
    # Initial treatment value derived from probability, p.s for antibiotic.s 
    
    
    for (i in 2:n.day){
        for (j in 1:n.bed){
            if(is.na(matrix_AtbTrt[i,j]) & (matrix_AtbTrt[i-1,j] != 0)){
                matrix_AtbTrt[i, j] <-  matrix_AtbTrt[i-1,j] + 1
            }else if (!is.na(matrix_AtbTrt[i,j])){
                matrix_AtbTrt[i, j] <-  matrix_AtbTrt[i, j]
            }else{
                matrix_AtbTrt[i, j] <-  0
            }
        }
    }
    #Case1 (patients with missing value and having antibiotic the day before): add one more day of days in the hospital
    #Case2 (patients with no missing value): the same value as it is
    #Case3 (all other patients i.e. without antibiotic or value=0): 0
    #Complete the matrix of cummulative length of stay for treated patients
    
    matrix_AtbTrt2 <- matrix(NA, nrow=n.day, ncol=n.bed)
    for (i in 1:n.day){
        for (j in 1:n.bed){
            if(matrix_DuraDay[i,j]>=matrix_AtbTrt[i,j] & matrix_AtbTrt[i,j]!=0){
                matrix_AtbTrt2[i, j] <-  1
            }else{
                matrix_AtbTrt2[i, j] <-  0
            }
        }
    }
    #Case1 (treated patients with atb duration longer than days staying in the hospital so far): 1
    #Case2 (untreated patients or duration of atb shorter than current days in the hospital so far): 0
    #Output matrix matrix_AtbTrt2 containing binary variable (treated vs not treated) for any bed on any particular day
    
    ## antibiotics to treat resistant organisms 
    matrix_DuraDay2 <- matrix(NA, nrow=n.day, ncol=n.bed)
    #matrix_DuraDay is the matrix with number of days of antibiotics.r for every patient in patient.matrix
    
    for (i in 1:max(patient.matrix)){
        for (j in 1:n.bed){
            matrix_DuraDay2[get1stdayofstay(i,j,patient.matrix), j]<-abs(round(rnorm(1, mean=meanDur.r*timestep, sd=sdDur*timestep)))
        }
    }
    #number of days of antibiotic.r is randomly drawn from a normal dist
    
    for (i in 2:n.day){
        for (j in 1:n.bed){
            if(is.na(matrix_DuraDay2[i,j])){
                matrix_DuraDay2[i, j] <- matrix_DuraDay2[i-1,j]
            }else{
                matrix_DuraDay2[i, j] <- matrix_DuraDay2[i, j]
            }
        }
    }
    #Fill the matrix for antibiotics.r with same dimension as patient.matrix
    
    matrix_AtbTrt.r <- matrix(NA, nrow=n.day, ncol=n.bed) 
    # #matrix_AtbTrt.r to count the cummulative length of stay for treated patients in the same dimension as patient.matrix
    
    for (i in 1:max(patient.matrix)){
        for (j in 1:n.bed){
            rand <- runif(1,0,1)
            if (rand < p.r.day1) {
                matrix_AtbTrt.r[get1stdayofstay(i,j,patient.matrix), j] <-  2
            } else {
                matrix_AtbTrt.r[get1stdayofstay(i,j,patient.matrix), j] <-  0
            }
        }
    }
    # Initial treatment value derived from probability, p.r.day1 for antibiotic.r 
    
    howmanydaysofabt<- function(m, i, j){ 
        n <- 0
        id <- patient.matrix[i, j]
        for (q in 1:(i-1)) {
            if (m[(i-q), j] > 0) {
                if (patient.matrix[(i-q), j] == id){
                    n <- n+1
                }else{
                    break
                }
            }
        }
        return (n)
    }

    for (i in 2:nrow(patient.matrix)){
        for (j in 1:n.bed){
            #print(paste("i", i, "j", j))
            #print(matrix_AtbTrt.r[i, j])
            if(is.na(matrix_AtbTrt.r[i, j])){
                rand <- runif(1,0,1)
                # case of no antibiotics for resistant organisms in the day before
                if (matrix_AtbTrt.r[i-1, j] == 0) {
                    if (rand < p.r.dayafter) {
                        matrix_AtbTrt.r [i,j] <- 2
                    } else {
                        matrix_AtbTrt.r [i,j] <- 0
                    }
                    # case of there is antibiotics for resistant organisms in the day before
                } else if (matrix_AtbTrt.r[i-1, j] == 2) { 
                    noabtday <- matrix_DuraDay2[i-1, j]
                    if (howmanydaysofabt(matrix_AtbTrt.r, i, j) < noabtday) {
                        matrix_AtbTrt.r [i,j] <- 2
                    } else {
                        matrix_AtbTrt.r [i,j] <- 0
                    }
                } else {
                    matrix_AtbTrt.r [i,j] <- 99 # Error code
                }
            }
        }
    }
    
    matrix_AtbTrt3<- matrix_AtbTrt2+matrix_AtbTrt.r
    #print('gen abx complete')
    
    return(matrix_AtbTrt3)
}

#################3. Generate baseline carriage status 

colo.table <- function(patient.matrix, los.array, prob_StartBact_R, prop_S_nonR, prop_Sr_inR, prop_sr_inR){
    
    prob_start_S <- prop_S_nonR*(1-prob_StartBact_R)
    prob_start_ss <- 1-prob_start_S-prob_StartBact_R
    prob_start_Sr <- prop_Sr_inR*prob_StartBact_R
    prob_start_sr <- prop_sr_inR*prob_StartBact_R
    prob_start_sR <- prob_StartBact_R-prob_start_Sr-prob_start_sr
    prob_StartBact_bi <- c(prob_start_S,prob_start_Sr,prob_start_sR,prob_start_sr)

    #Generating a vector of random status with runif (change for other distribution)
    number_of_patients <- dim(los.array)[2]
    Patient_unif <- runif(number_of_patients,0,1)
    Patient_StartBact <- rep(NA, number_of_patients)
    Patient_StartBact[Patient_unif > sum(prob_StartBact_bi)] <- 'ss'
    Patient_StartBact[(Patient_unif > sum(prob_StartBact_bi[1:3])) & (Patient_unif <= sum(prob_StartBact_bi))] <- 'sr'
    Patient_StartBact[(Patient_unif > sum(prob_StartBact_bi[1:2])) & (Patient_unif <= sum(prob_StartBact_bi[1:3]))] <- 'sR'
    Patient_StartBact[(Patient_unif > sum(prob_StartBact_bi[1])) & (Patient_unif <= sum(prob_StartBact_bi[1:2]))] <- 'Sr'
    Patient_StartBact[Patient_unif <= sum(prob_StartBact_bi[1])] <- 'S'
    
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

####################4. Update values for every day  
nextDay <- function(patient.matrix, abx.matrix, colo.matrix, 
                    pi_r1, bif, mu1, mu2, repop.r1, repop.r2, 
                    repop.s1, repop.s2,repop.s3, abx.r, abx.s){
    
    pi_r2 <- pi_r1 * bif                 # pi_r2= probability of R transmitting to s to become sr 
    
    # For each day (first day should be filled)
    for(i in 2:nrow(patient.matrix)){
        # For each bed
        for(j in 1:ncol(patient.matrix)){
            #case S
            #print(paste("i:", i, "j:", j))
            #print(colo.matrix[i-1, j])
            if(is.na(colo.matrix[i, j])){
                r_num <- sum(colo.matrix[i-1,] == "sR") #only R can be transmitted 
                if(colo.matrix[i-1, j] == "S"){
                    #print("----case S")
                    # check antibiotic
                    roll_clear <- runif(1, 0, 1)
                    roll_transmit <- runif(1, 0, 1)
                    prob_r <- 1-((1-pi_r1)^r_num)
                    if (abx.matrix[i-1, j] == 1 & roll_clear < abx.s){
                        colo.matrix[i, j] <- "ss"
                    } else if (abx.matrix[i-1, j] > 1 & roll_clear < abx.r){
                        colo.matrix[i, j] <- "ss"
                    } else if (roll_transmit < prob_r){ 
                        colo.matrix[i, j] <- "Sr"
                    }else {
                        colo.matrix[i, j] <- "S"
                    }
                    
                    # case ss
                }else if(colo.matrix[i-1, j] == "ss"){
                    #print("----case s")
                    # roll for transmission of r
                    roll_r <- runif(1, 0, 1)
                    prob_r <- 1-((1-pi_r2)^r_num)
                    # roll for repopulation of s to become S
                    roll_ss <- runif(1, 0, 1)
                    if (roll_r < prob_r) { 
                        colo.matrix[i,j]<-"sr" 
                    } else if ( roll_ss < repop.s1) {
                        colo.matrix[i,j]<-"S" 
                    } else{
                        colo.matrix[i, j] <- "ss"
                    }
                    
                    # case Sr
                }else if(colo.matrix[i-1, j] == "Sr"){
                    #print("----case Sr")
                    # check antibiotics 
                    roll_clear <- runif(1, 0, 1)
                    roll_decolonise <- runif(1, 0, 1)
                    if(abx.matrix[i-1, j] ==1 & roll_clear < abx.s){
                        colo.matrix[i, j] <- "sr"
                    } else if (abx.matrix[i-1, j] >1 & roll_clear < abx.r ) {
                        colo.matrix[i, j] <- "sr"
                    } else if(roll_decolonise < mu1){ 
                        colo.matrix[i, j] <- "S"
                    }else {
                        colo.matrix[i, j] <- "Sr"
                    }
                    
                    # case sr
                }else if(colo.matrix[i-1, j] == "sr"){
                    #print("----case sr")
                    roll_repop <- runif(1, 0, 1)
                    roll_decolonise <- runif(1, 0, 1)
                    # check antibiotics
                    if(abx.matrix[i-1, j] == 1 & roll_repop < repop.r2){
                        colo.matrix[i, j] <- "sR"
                    }else if(abx.matrix[i-1, j] == 0 & roll_repop < repop.s2){
                        colo.matrix[i, j] <- "Sr"
                        # }else if(abx.matrix[i-1, j] == 0 & roll_repop < repop.r3){
                        #     colo.matrix[i, j] <- "sR"
                    }else if (roll_decolonise < mu2){
                        colo.matrix[i, j] <- "ss"
                    }else {
                        colo.matrix[i, j] <- "sr"
                    }
                    
                    # case sR
                }else if(colo.matrix[i-1, j] == "sR"){
                    #print("----case sR")
                    roll_clear <- runif(1, 0, 1)
                    if(abx.matrix[i-1, j] > 1 & roll_clear < abx.r){
                        colo.matrix[i, j] <- "sr"
                    }else if (abx.matrix[i-1, j] == 0 & roll_clear < repop.s3){
                        colo.matrix[i, j] <- "sr"
                    }else {
                        colo.matrix[i, j] <- "sR"
                    }
                }else{
                    print("error")
                    colo.matrix[i, j] <- "E"
                }
            } # if 0
        }  # for j
    } # for i
    
    return(colo.matrix)
}

diff_prevalence <- function(n.bed, mean.max.los, p.s, p.r.day1, p.r.dayafter,
                            prob_StartBact_R, prop_S_nonR, prop_Sr_inR, prop_sr_inR,
                            pi_r1, bif, mu1, mu2, abx.r, abx.s,
                            repop.r1, repop.r2, repop.s1, repop.s2, repop.s3,
                            short_dur, long_dur){
    n.day <- 30
    iterations <- 10
    iter_totalsR <- matrix(NA, nrow = n.day, ncol = iterations)
    
    for(iter in 1:iterations){
        
        #print(paste("iter:", iter, "y:", y_count, '-', y, "x", x_count, '-', x))
        #Generate length of stay and antibiotic duration table
        abx_iter <- abx.table(n.bed=n.bed, n.day=n.day, mean.max.los=mean.max.los, p.s=p.s, p.r.day1=p.r.day1, p.r.dayafter=p.r.dayafter, meanDur=short_dur)
        #Generate baseline carriage status
        array_LOS_iter <- array_LOS_func(los_duration=abx_iter[[2]])
        #Update values for every day
        array_StartBact_iter <- gen_StartBact(los=array_LOS_iter, prob_StartBact_R=prob_StartBact_R, 
                                              prop_S_nonR=prop_S_nonR, prop_Sr_inR=prop_Sr_inR, 
                                              prop_sr_inR=prop_sr_inR, 
                                              n.bed=n.bed, n.day=n.day)
        #output
        colo.matrix_filled_iter <- nextDay(patient.matrix= abx_iter[[1]], array_LOS=array_LOS_iter, 
                                          abx.matrix=abx_iter[[3]], colo.matrix=array_StartBact_iter, 
                                          pi_r1=pi_r1, bif=bif, mu1=mu1, mu2=mu2, 
                                          abx.r=abx.r,abx.s=abx.s,
                                          repop.r1 = repop.r1, repop.r2 = repop.r2,
                                          repop.s1 = repop.s1, repop.s2 = repop.s2,repop.s3 = repop.s3)
        #Summary
        df <- data.frame(colo.matrix_filled_iter)
        iter_totalsR[, iter] <- rowSums(df == "sR")
        #print("end iteration loop")
    }
    totalsR_short <- mean(rowSums(iter_totalsR[ceiling(n.day*1/3):nrow(iter_totalsR),])/iterations/n.bed)
    
    iter_totalsR <- matrix(NA, nrow = n.day, ncol = iterations)
    for(iter in 1:iterations){
        
        #print(paste("iter:", iter, "y:", y_count, '-', y, "x", x_count, '-', x))
        #Generate length of stay and antibiotic duration table
        abx_iter <- abx.table(n.bed=n.bed, n.day=n.day, mean.max.los=mean.max.los, p.s=p.s, p.r.day1=p.r.day1, p.r.dayafter=p.r.dayafter, meanDur=long_dur)
        #Generate baseline carriage status
        array_LOS_iter <- array_LOS_func(los_duration=abx_iter[[2]])
        #Update values for every day
        array_StartBact_iter <- gen_StartBact(los=array_LOS_iter,prob_StartBact_R=prob_StartBact_R, 
                                              prop_S_nonR=prop_S_nonR, prop_sr_inR=prop_sr_inR,
                                              prop_Sr_inR=prop_Sr_inR, n.bed=n.bed, n.day=n.day)
        #output
        colo.matrix_filled_iter <- nextDay(patient.matrix= abx_iter[[1]], array_LOS=array_LOS_iter, 
                                          abx.matrix=abx_iter[[3]], colo.matrix=array_StartBact_iter, 
                                          pi_r1=pi_r1, bif=bif, mu1=mu1, mu2=mu2, 
                                          abx.r=abx.r,abx.s=abx.s,
                                          repop.r1 = repop.r1, repop.r2 = repop.r2,
                                          repop.s1 = repop.s1, repop.s2 = repop.s2, repop.s3 = repop.s3)
        #Summary
        df <- data.frame(colo.matrix_filled_iter)
        iter_totalsR[,iter] <- rowSums(df == "sR")
        #print("end iteration loop")
    }
    totalsR_long <- mean(rowSums(iter_totalsR[ceiling(n.day*1/3):nrow(iter_totalsR),])/iterations/n.bed)
    
    #print(paste("totalsR_long", totalsR_long, "totalsR_short", totalsR_short))
    
    return(totalsR_long - totalsR_short)
}

# diff_prevalence(n.bed=20, mean.max.los=4, p.s=0.1, p.r.day1=0.2, p.r.dayafter=0.01,
#                 prob_StartBact_R=0.3, prop_S_nonR=0.1, prop_Sr_inR=0.1, prop_sr_inR=0.1, 
#                 pi_r1=0.1, bif=2, mu1=.1, mu2=.1, abx.r=.1, abx.s=.1,
#                 repop.r1=.1, repop.r2=.1, repop.s1=.1, repop.s2=.1, repop.s3=.1,
#                 short_dur=2, long_dur=10)

