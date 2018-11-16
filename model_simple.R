# positive_norm_sample <- function(mean, sd){                      
#     v<-round(rnorm(1, mean=mean, sd=sd))
#     if (v < 0) {
#         v<-v*-1
#     }
#     return(v)
# }

abx.table<- function(n.bed, n.days, mean.max.los, p, meanDur) {
    
    # generate a table of number of days we want to observe (rows) -
    # against number of beds in the ward (columns), filled in with patient id numbers
    
    patient.matrix <- matrix(NA, nrow=n.days, ncol=n.bed)
    #make up a matrix of number of days we want to observe (rows) -
    #against number of beds in the ward (columns)
    
    n.patient <- n.bed*n.days 
    #generate patient id numbers, the maximum number of patients possible is number of bed multiple by
    #number of days. This is to ensure there are enough total number of patients generated to fill table
    
    patient.id <- c(1:n.patient)    
    #vectorise the patient id to be used for filling in the patient.matrix
    
    vector.los <- rep(0,n.days)
    #Generating a vector with 0s 
    #make up an empty vector for 0th column of matrix so that the next column starts with 
    #(last patient number of the previous column+i)
    
    for (j in 1:n.bed) {
        #for each column (representing different beds in the ward) in patient.matrix
        
        los <- c()
        final <- vector.los[n.days]
        #last patient number in previous column
        
        for (i in 1:n.days) {
            #for each row (representing number of days we want to observe)
            
            los <- rep(patient.id[i+final], rexp(1,1/mean.max.los))
            #repeat the (last patient number of the previous column+i) random number of times
            
            vector.los <- c(vector.los, los)
            #combine 0th column (los0) and subsequent columns (los) into one vector
        }
        
        vector.los<- vector.los[-(1:n.days)]
        #remove 0th column from patient.matrix
        
        if (length(vector.los) > n.days) {                     
            vector.los<-vector.los[1:n.days]
        }
        #ensure each column is of the same length as number of days we are observing 
        #(number of rows)
        
        patient.matrix[,j] <- vector.los
        #fill the columns (different beds in the ward) with length of stay vector generated
    }
    
    los_duration<-as.vector(patient.matrix)
    #frequency summary of patient.matrix - patient id against number of days of stay for each patient
    
    #generate antibiotic use table
    
    # matrix_DuraDay is the matrix with number of days of antibiotics for every patient in patient.matrix
    matrix_DuraDay <- matrix(NA, nrow=n.days, ncol=n.bed)
    #returns first day of stay for patientnumber in bednumber, or NA if patient not there
    get1stdayofstay<-function(patientnumber,bednumber, bedoccmat){
        match(patientnumber, bedoccmat[,bednumber])
    }
    
    #number of days of antibiotic is randomly drawn from a uniform dist with min=minDur, max=maxDur
    for (i in 1:max(patient.matrix)){
        for (j in 1:n.bed){
            matrix_DuraDay[get1stdayofstay(i,j,patient.matrix), j]<-abs(round(rnorm(1, mean=meanDur, sd=5)))
        }
    }
    
    # Fill the matrix with same dimension as patient.matrix
    for (i in 2:n.days){
        for (j in 1:n.bed){
            matrix_DuraDay[i, j] <-  if(is.na(matrix_DuraDay[i,j])){
                matrix_DuraDay[i-1,j]
            }
            else{matrix_DuraDay[i, j]}
        }
    }
    
    # matrix_AtbTrt to count the cummulative length of stay for treated patients in the same dimension as patient.matrix
    matrix_AtbTrt <- matrix(NA, nrow=n.days, ncol=n.bed) # empty matrix
    # Initial treatment value derived from probability, p
    for (i in 1:max(patient.matrix)){
        for (j in 1:n.bed){
            rand <- runif(1,0,1)
            if (rand < p) {
                matrix_AtbTrt[get1stdayofstay(i,j,patient.matrix), j] <-  1
            } else {
                matrix_AtbTrt[get1stdayofstay(i,j,patient.matrix), j] <-  0
            }
        }
    }
    # Initial treatment value derived from probability, p for narrow spectrum antibiotic
    
    for (i in 2:n.days){
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
    
    # Output matrix matrix_AtbTrt2 containing binary variable (treated vs not treated) for any bed on any particular day
    matrix_AtbTrt2 <- matrix(NA, nrow=n.days, ncol=n.bed)
    for (i in 1:n.days){
        for (j in 1:n.bed){
            if(matrix_DuraDay[i,j]>=matrix_AtbTrt[i,j] & matrix_AtbTrt[i,j]!=0){
                matrix_AtbTrt2[i, j] <-  1
            }else{
                matrix_AtbTrt2[i, j] <-  0
            }
        }
    }
    
    return(list(patient.matrix, los_duration, matrix_AtbTrt2))
}

array_LOS_func<- function(los_duration) {
    los.dur<-as.vector(table(los_duration))
    array_LOS<-array(dim=c(2,length(los.dur)))
    array_LOS[1,]<-c(1:length(los.dur))
    array_LOS[2,]<-los.dur
    
    return(array_LOS)
}

#################3. Generate baseline carriage status ##################

gen_StartBact <- function(los, prob_StartBact_R, prop_S_nonR, n.bed, n.days){
    
    #define probabilities of importing Sensitive(S) or Resistant(R) bacteria, or nothing (N)
    prob_start_S <- prop_S_nonR*(1-prob_StartBact_R)
    prob_StartBact <- c(prob_start_S,prob_StartBact_R)
    
    stopifnot(sum(prob_StartBact) < 1) # Assert all probabilities combined are less than 1
    
    #Generating a vector of random status with runif (change for other distribution)
    number_of_patients <- dim(los)[2]
    Patient_unif <- runif(number_of_patients,0,1)
    Patient_StartBact <- rep(NA, number_of_patients)
    Patient_StartBact[Patient_unif > (prob_start_S+prob_StartBact_R)] <- 'N'
    Patient_StartBact[(Patient_unif <= prob_start_S+prob_StartBact_R) & (Patient_unif > prob_start_S)] <- 'R'
    Patient_StartBact[Patient_unif <= prob_start_S] <- 'S'
    
    #Creating array for carriage status
    array_StartBact <- matrix(NA, n.days, n.bed)
    
    end_idx <- 1
    for(i in 1:number_of_patients){
        # print(i)
        # print(paste(end_idx, end_idx + los[2, i] - 1))
        # print(c(Patient_StartBact[i], rep(NA, los[2, i]-1)))
        # print(length(c(Patient_StartBact[i], rep(NA, los[2, i]-1))))
        array_StartBact[end_idx:(end_idx + los[2, i] - 1)] <- c(Patient_StartBact[i], rep(NA, los[2, i]-1))
        end_idx = end_idx + los[2, i]
    }
    
    return(array_StartBact)
}

####################4. Update values for every day  #####################
nextDay <- function(bed_table, array_LOS, treat_table, colo_table, pi_sr, mu_r, pi_s, pi_r, abx.clear){
    
    # For each day (first day should be filled)
    for(i in 2:nrow(bed_table)){
        # For each bed
        for(j in 1:ncol(bed_table)){
            
            # case R
            #print(paste("i:", i, "j:", j))
            #print(colo_table[i-1, j])
            if(is.na(colo_table[i, j])){
                if(colo_table[i-1, j] == "R"){
                    #print("----case R")
                    # roll for clearance
                    roll <- runif(1, 0, 1)
                    if(roll < mu_r){
                        colo_table[i, j] <- "S"
                    }else{
                        colo_table[i, j] <- "R"
                    }
                    # case S
                }else if(colo_table[i-1, j] == "S"){
                    #print("----case S")
                    # check antibiotic
                    roll_clear <- runif(1, 0, 1)
                    roll_select <- runif(1, 0, 1)
                    if(treat_table[i-1, j] > 0 & roll_clear < abx.clear){
                        colo_table[i, j] <- "N"
                    # }else if(roll_clear < mu_s){
                    #     # roll for clearance
                    #     colo_table[i, j] <- "N"
                    }else if(roll_select < pi_sr){
                        # roll for selection
                        colo_table[i, j] <- "R"
                    }else{
                        colo_table[i, j] <- "S"
                    }
                    # case N
                }else if(colo_table[i-1, j] == "N"){
                    #print("----case N")
                    # roll for transmission of R
                    roll_r <- runif(1, 0, 1)
                    r_num <- sum(colo_table[i-1,] == "R")
                    prob_r <- 1-((1-pi_r)^r_num)
                    transmit_r <- roll_r < prob_r
                    # roll for transmission of S
                    roll_s <- runif(1, 0, 1)
                    s_num <- sum(colo_table[i-1,] == "S")
                    prob_s <- 1-((1-pi_s)^s_num)
                    transmit_s <- roll_s < prob_s
                    
                    if(transmit_r & transmit_s){
                        # if both R and S pass, randomly pick one of them
                        roll <- runif(1, 0, 1)
                        if(roll > 0.5| treat_table[i-1,j]){
                            colo_table[i, j] <- "R"
                        }else{
                            colo_table[i,j]<-"S"}
                    }else if(transmit_r & !transmit_s){
                        colo_table[i, j] <- "R"
                    }else if(!transmit_r & transmit_s& !treat_table[i-1,j]){
                        colo_table[i, j] <- "S"
                    }else{
                        colo_table[i, j] <- "N"
                    }
                }else{
                    print("error")
                    colo_table[i, j] <- "E"
                }
            }
        }  # for j
    } # for i
    
    return(colo_table)
}

diff_prevalence <- function(n.bed, mean.max.los, p, prop_S_nonR, 
                            prob_StartBact_R, pi_sr, mu_s, mu_r, pi_s, pi_r, abx.clear,
                            short_dur, long_dur){
    n.days<-30
    iterations<-10
    iter_totalR <- matrix(NA, nrow = n.days, ncol = iterations)
    for(iter in 1:iterations){
        
        abx_iter<-abx.table(n.bed, n.days, mean.max.los, p, meanDur = short_dur)
        
        array_LOS_iter<-array_LOS_func(los_duration=abx_iter[2])
        
        array_StartBact_iter<-gen_StartBact(los=array_LOS_iter, prob_StartBact_R=prob_StartBact_R,prop_S_nonR=prop_S_nonR,n.days=n.days, n.bed=n.bed)
        
        colo_table_filled_iter <- nextDay(bed_table= abx_iter[[1]], array_LOS=array_LOS_iter, treat_table=abx_iter[[3]], 
                                           colo_table=array_StartBact_iter, pi_sr=pi_sr, mu_s=mu_s, mu_r=mu_r, pi_s=pi_s, pi_r=pi_r, abx.clear=abx.clear)
    
        #Summary
        df <- data.frame(colo_table_filled_iter)
        iter_totalR[, iter] <- rowSums(df == "R")    
    }
    totalR_short <- mean(rowSums(iter_totalR)/iterations/n.bed)
    
    iter_totalR <- matrix(NA, nrow = n.days, ncol = iterations)
    for(iter in 1:iterations){
        
        abx_iter<-abx.table(n.bed, n.days, mean.max.los, p, meanDur = long_dur)
        
        array_LOS_iter<-array_LOS_func(los_duration=abx_iter[2])
        
        array_StartBact_iter<-gen_StartBact(los=array_LOS_iter, prob_StartBact_R=prob_StartBact_R,prop_S_nonR=prop_S_nonR, n.days=n.days, n.bed=n.bed)
        
        colo_table_filled_iter <- nextDay(bed_table= abx_iter[[1]], array_LOS=array_LOS_iter, treat_table=abx_iter[[3]], 
                                          colo_table=array_StartBact_iter, pi_sr=pi_sr, mu_s=mu_s, mu_r=mu_r, pi_s=pi_s, pi_r=pi_r, abx.clear=abx.clear)
        
        #Summary
        df <- data.frame(colo_table_filled_iter)
        iter_totalR[, iter] <- rowSums(df == "R")    
    }
    totalR_long <- mean(rowSums(iter_totalR)/iterations/n.bed)
    
    print(paste("totalR_long", totalR_long, "totalR_short", totalR_short))
    
    return(totalR_long - totalR_short)
}

diff_prevalence(n.bed =20, mean.max.los=3, p=0.1,prob_StartBact_R=0.4, prop_S_nonR=0.3,
                pi_sr=0.1, mu_s=0.2, mu_r=0.1, pi_s=0.1, pi_r=0.1, abx.clear=0.1, short_dur = 4, long_dur = 10)

