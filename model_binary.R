# ######## Functions ########
# positive_norm_sample <- function(mean, sd){                      
#     v<-round(rnorm(1, mean=mean, sd=sd))
#     if (v < 0) {
#         v<-v*-1
#     }
#     return(v)
# }
#########2. Generate length of stay and antibiotic duration table
#(1 for short duration and 1 for long duration) 
#simulate inpatients with various lengths of stay
#allocate various duration of antibiotics for each patient

abx.table<- function (n.bed, n.days, mean.max.los, p.s, p.r.day1, p.r.dayafter, meanDur) {
    # generate a table of number of days we want to observe (rows) -
    # against number of beds in the ward (columns), filled in with patient id numbers
    
    #make up a matrix of number of days we want to observe (rows) -
    #against number of beds in the ward (columns)
    patient.matrix <- matrix(NA, nrow=n.days, ncol=n.bed)
   
    #generate patient id numbers, the maximum number of patients possible is number of bed multiple by
    #number of days. This is to ensure there are enough total number of patients generated to fill table
    n.patient <- n.bed*n.days 
    
    #vectorise the patient id to be used for filling in the patient.matrix
    patient.id <- c(1:n.patient)    
    
    #Generating a vector with 0s 
    #make up an empty vector for 0th column of matrix so that the next column starts with 
    #(last patient number of the previous column+i)
    vector.los <- rep(0,n.days)

    #for each column (representing different beds in the ward) in patient.matrix
    for (j in 1:n.bed) {
        
        los <- c()
        final <- vector.los[n.days]
        #last patient number in previous column
        
        for (i in 1:n.days) {
            #for each row (representing number of days we want to observe)
            
            los <- rep(patient.id[i+final], as.integer(rexp(1, 1/mean.max.los)))
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
    ## antibiotics to treat sensitive organisms 
    matrix_DuraDay <- matrix(NA, nrow=n.days, ncol=n.bed)
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
            matrix_DuraDay[get1stdayofstay(i,j,patient.matrix), j] <- abs(round(rnorm(1, mean=meanDur, sd=5)))
        }
    }
    #number of days of antibiotic.s is randomly drawn from a normal dist
    
    for (i in 2:n.days){
        for (j in 1:n.bed){
            if(is.na(matrix_DuraDay[i,j])){
                matrix_DuraDay[i, j] <- matrix_DuraDay[i-1,j]
            }else{
                matrix_DuraDay[i, j] <- matrix_DuraDay[i, j]
            }
        }
    }
    #Fill the matrix for antibiotics.s with same dimension as patient.matrix
    
    matrix_AtbTrt <- matrix(NA, nrow=n.days, ncol=n.bed) 
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
    #Case1 (patients with missing value and having antibiotic the day before): add one more day of days in the hospital
    #Case2 (patients with no missing value): the same value as it is
    #Case3 (all other patients i.e. without antibiotic or value=0): 0
    #Complete the matrix of cummulative length of stay for treated patients
    
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
    #Case1 (treated patients with atb duration longer than days staying in the hospital so far): 1
    #Case2 (untreated patients or duration of atb shorter than current days in the hospital so far): 0
    #Output matrix matrix_AtbTrt2 containing binary variable (treated vs not treated) for any bed on any particular day
    
    ## antibiotics to treat resistant organisms 
    matrix_DuraDay2 <- matrix(NA, nrow=n.days, ncol=n.bed)
    #matrix_DuraDay is the matrix with number of days of antibiotics.r for every patient in patient.matrix
    
    for (i in 1:max(patient.matrix)){
        for (j in 1:n.bed){
            matrix_DuraDay2[get1stdayofstay(i,j,patient.matrix), j]<-abs(round(rnorm(1, mean=meanDur, sd=5)))
        }
    }
    #number of days of antibiotic.r is randomly drawn from a normal dist
    
    for (i in 2:n.days){
        for (j in 1:n.bed){
            if(is.na(matrix_DuraDay2[i,j])){
                matrix_DuraDay2[i, j] <- matrix_DuraDay2[i-1,j]
            }else{
                matrix_DuraDay2[i, j] <- matrix_DuraDay2[i, j]
            }
        }
    }
    #Fill the matrix for antibiotics.r with same dimension as patient.matrix
    
    matrix_AtbTrt.r <- matrix(NA, nrow=n.days, ncol=n.bed) 
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
    # function to check how many days the patient has been on antibiotics 
    
    # Unit test for function howmanydaysofabt
    # patient.matrix <- c(1, 1, 2, 2, 3, 3, 4, 4, 4, 5, 5, 5, 6, 7, 7, 7, 7, 7)
    # dim(patient.matrix) <- c(6, 3)
    # m <- c(0, 2, 2, 2, 2, 0, 2, 0, 0, 2, 2, 0, 2, 2, 2, 2, 2, 0)
    # dim(m) <- c(6, 3)
    # Check base case
    # howmanydaysofabt(m, 2, 1) should return 2
    # howmanydaysofabt(m, 2, 2) should return 2
    # Check cross over between patients
    # howmanydaysofabt(m, 4, 1) should return 3
    # howmanydaysofabt(m, 4, 3) should return 4
    
    
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
    
    return(list(patient.matrix, los_duration, matrix_AtbTrt3))
}

array_LOS_func<- function(los_duration) {
    los.dur<-as.vector(table(los_duration))
    array_LOS<-array(dim=c(2,length(los.dur)))
    array_LOS[1,]<-c(1:length(los.dur))
    array_LOS[2,]<-los.dur
    
    return(array_LOS)
}


#################3. Generate baseline carriage status 

gen_StartBact <- function(los, prob_StartBact_R, prop_S_nonR, prop_Sr_inR, prop_sr_inR, n.bed, n.days){
    
    prob_start_S <- prop_S_nonR*(1-prob_StartBact_R)
    prob_start_s <- 1-prob_start_S-prob_StartBact_R
    prob_start_Sr <- prop_Sr_inR*prob_StartBact_R
    prob_start_sr <- prop_sr_inR* prob_StartBact_R
    prob_start_sR <- (1-prop_Sr_inR-prop_sr_inR)*prob_StartBact_R
    prob_StartBact_bi <- c(prob_start_S,prob_start_Sr,prob_start_sR,prob_start_sr)
    
    stopifnot(sum(prob_StartBact_bi) < 1) # Assert all probabilities combined are less than 1

    #Generating a vector of random status with runif (change for other distribution)
    number_of_patients <- dim(los)[2]
    Patient_unif <- runif(number_of_patients, 0, 1)
    Patient_StartBact <- rep(NA, number_of_patients)
    Patient_StartBact[Patient_unif > sum(prob_StartBact_bi)] <- 's'
    Patient_StartBact[(Patient_unif > sum(prob_StartBact_bi[1:3])) & (Patient_unif <= sum(prob_StartBact_bi))] <- 'sr'
    Patient_StartBact[(Patient_unif > sum(prob_StartBact_bi[1:2])) & (Patient_unif <= sum(prob_StartBact_bi[1:3]))] <- 'sR'
    Patient_StartBact[(Patient_unif > sum(prob_StartBact_bi[1])) & (Patient_unif <= sum(prob_StartBact_bi[1:2]))] <- 'Sr'
    Patient_StartBact[Patient_unif <= sum(prob_StartBact_bi[1])] <- 'S'
    
    #Creating matrix for carriage status
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
    
    #print('gen bact complete')
    return(array_StartBact)
}



####################4. Update values for every day  
nextDay <- function(bed_table, array_LOS, treat_table, colo_table, 
                    pi_r1, bif, mu1, mu2, repop.r1, repop.r2, 
                    repop.s1, repop.s2,repop.s3, abx.r, abx.s){
    
    pi_r2 <- pi_r1 * bif                 # pi_r2= probability of R transmitting to s to become sr 
    
    # For each day (first day should be filled)
    for(i in 2:nrow(bed_table)){
        # For each bed
        for(j in 1:ncol(bed_table)){
            #case S
            #print(paste("i:", i, "j:", j))
            #print(colo_table[i-1, j])
            if(is.na(colo_table[i, j])){
                r_num <- sum(colo_table[i-1,] == "sR") #only R can be transmitted 
                if(colo_table[i-1, j] == "S"){
                    #print("----case S")
                    # check antibiotic
                    roll_clear <- runif(1, 0, 1)
                    roll_transmit <- runif(1, 0, 1)
                    prob_r <- 1-((1-pi_r1)^r_num)
                    if (treat_table[i-1, j] == 1 & roll_clear < abx.s){
                        colo_table[i, j] <- "s"
                    } else if (treat_table[i-1, j] > 1 & roll_clear < abx.r){
                        colo_table[i, j] <- "s"
                    } else if (roll_transmit < prob_r){ 
                        colo_table[i, j] <- "Sr"
                    }else {
                        colo_table[i, j] <- "S"
                    }
                    
                    # case s
                }else if(colo_table[i-1, j] == "s"){
                    #print("----case s")
                    # roll for transmission of r
                    roll_r <- runif(1, 0, 1)
                    prob_r <- 1-((1-pi_r2)^r_num)
                    # roll for repopulation of s to become S
                    roll_s <- runif(1, 0, 1)
                    if (roll_r < prob_r) { 
                        colo_table[i,j]<-"sr" 
                    } else if ( roll_s < repop.s1) {
                        colo_table[i,j]<-"S" 
                    } else{
                        colo_table[i, j] <- "s"
                    }
                    
                    # case Sr
                }else if(colo_table[i-1, j] == "Sr"){
                    #print("----case Sr")
                    # check antibiotics 
                    roll_clear <- runif(1, 0, 1)
                    roll_decolonise <- runif(1, 0, 1)
                    if(treat_table[i-1, j] ==1 & roll_clear < abx.s){
                        colo_table[i, j] <- "sr"
                    } else if (treat_table[i-1, j] >1 & roll_clear < abx.r ) {
                        colo_table[i, j] <- "sr"
                    } else if(roll_decolonise < mu1){ 
                        colo_table[i, j] <- "S"
                    }else {
                        colo_table[i, j] <- "Sr"
                    }
                    
                    # case sr
                }else if(colo_table[i-1, j] == "sr"){
                    #print("----case sr")
                    roll_repop <- runif(1, 0, 1)
                    roll_decolonise <- runif(1, 0, 1)
                    # check antibiotics
                    if(treat_table[i-1, j] == 1 & roll_repop < repop.r2){
                        colo_table[i, j] <- "sR"
                    }else if(treat_table[i-1, j] == 0 & roll_repop < repop.s2){
                        colo_table[i, j] <- "Sr"
                        # }else if(treat_table[i-1, j] == 0 & roll_repop < repop.r3){
                        #     colo_table[i, j] <- "sR"
                    }else if (roll_decolonise < mu2){
                        colo_table[i, j] <- "s"
                    }else {
                        colo_table[i, j] <- "sr"
                    }
                    
                    # case sR
                }else if(colo_table[i-1, j] == "sR"){
                    #print("----case sR")
                    roll_clear <- runif(1, 0, 1)
                    if(treat_table[i-1, j] > 1 & roll_clear < abx.r){
                        colo_table[i, j] <- "sr"
                    }else if (treat_table[i-1, j] == 0 & roll_clear < repop.s3){
                        colo_table[i, j] <- "sr"
                    }else {
                        colo_table[i, j] <- "sR"
                    }
                }else{
                    print("error")
                    colo_table[i, j] <- "E"
                }
            } # if 0
        }  # for j
    } # for i
    
    return(colo_table)
}

diff_prevalence <- function(n.bed, mean.max.los, p.s, p.r.day1, p.r.dayafter,
                            prob_StartBact_R, prop_S_nonR, prop_Sr_inR, prop_sr_inR,
                            pi_r1, bif, mu1, mu2, abx.r, abx.s,
                            repop.r1, repop.r2, repop.s1, repop.s2, repop.s3,
                            short_dur, long_dur){
    n.days <- 30
    iterations <- 10
    iter_totalsR <- matrix(NA, nrow = n.days, ncol = iterations)
    
    for(iter in 1:iterations){
        
        #print(paste("iter:", iter, "y:", y_count, '-', y, "x", x_count, '-', x))
        #Generate length of stay and antibiotic duration table
        abx_iter <- abx.table(n.bed=n.bed, n.days=n.days, mean.max.los=mean.max.los, p.s=p.s, p.r.day1=p.r.day1, p.r.dayafter=p.r.dayafter, meanDur=short_dur)
        #Generate baseline carriage status
        array_LOS_iter <- array_LOS_func(los_duration=abx_iter[[2]])
        #Update values for every day
        array_StartBact_iter <- gen_StartBact(los=array_LOS_iter, prob_StartBact_R=prob_StartBact_R, 
                                              prop_S_nonR=prop_S_nonR, prop_Sr_inR=prop_Sr_inR, 
                                              prop_sr_inR=prop_sr_inR, 
                                              n.bed=n.bed, n.days=n.days)
        #output
        colo_table_filled_iter <- nextDay(bed_table= abx_iter[[1]], array_LOS=array_LOS_iter, 
                                          treat_table=abx_iter[[3]], colo_table=array_StartBact_iter, 
                                          pi_r1=pi_r1, bif=bif, mu1=mu1, mu2=mu2, 
                                          abx.r=abx.r,abx.s=abx.s,
                                          repop.r1 = repop.r1, repop.r2 = repop.r2,
                                          repop.s1 = repop.s1, repop.s2 = repop.s2,repop.s3 = repop.s3)
        #Summary
        df <- data.frame(colo_table_filled_iter)
        iter_totalsR[, iter] <- rowSums(df == "sR")
        #print("end iteration loop")
    }
    totalsR_short <- mean(rowSums(iter_totalsR[ceiling(n.days*1/3):nrow(iter_totalsR),])/iterations/n.bed)
    
    iter_totalsR <- matrix(NA, nrow = n.days, ncol = iterations)
    for(iter in 1:iterations){
        
        #print(paste("iter:", iter, "y:", y_count, '-', y, "x", x_count, '-', x))
        #Generate length of stay and antibiotic duration table
        abx_iter <- abx.table(n.bed=n.bed, n.days=n.days, mean.max.los=mean.max.los, p.s=p.s, p.r.day1=p.r.day1, p.r.dayafter=p.r.dayafter, meanDur=long_dur)
        #Generate baseline carriage status
        array_LOS_iter <- array_LOS_func(los_duration=abx_iter[[2]])
        #Update values for every day
        array_StartBact_iter <- gen_StartBact(los=array_LOS_iter,prob_StartBact_R=prob_StartBact_R, 
                                              prop_S_nonR=prop_S_nonR, prop_sr_inR=prop_sr_inR,
                                              prop_Sr_inR=prop_Sr_inR, n.bed=n.bed, n.days=n.days)
        #output
        colo_table_filled_iter <- nextDay(bed_table= abx_iter[[1]], array_LOS=array_LOS_iter, 
                                          treat_table=abx_iter[[3]], colo_table=array_StartBact_iter, 
                                          pi_r1=pi_r1, bif=bif, mu1=mu1, mu2=mu2, 
                                          abx.r=abx.r,abx.s=abx.s,
                                          repop.r1 = repop.r1, repop.r2 = repop.r2,
                                          repop.s1 = repop.s1, repop.s2 = repop.s2, repop.s3 = repop.s3)
        #Summary
        df <- data.frame(colo_table_filled_iter)
        iter_totalsR[,iter] <- rowSums(df == "sR")
        #print("end iteration loop")
    }
    totalsR_long <- mean(rowSums(iter_totalsR[ceiling(n.days*1/3):nrow(iter_totalsR),])/iterations/n.bed)
    
    #print(paste("totalsR_long", totalsR_long, "totalsR_short", totalsR_short))
    
    return(totalsR_long - totalsR_short)
}

# diff_prevalence(n.bed=20, mean.max.los=4, p.s=0.1, p.r.day1=0.2, p.r.dayafter=0.01,
#                 prob_StartBact_R=0.3, prop_S_nonR=0.1, prop_Sr_inR=0.1, prop_sr_inR=0.1, 
#                 pi_r1=0.1, bif=2, mu1=.1, mu2=.1, abx.r=.1, abx.s=.1,
#                 repop.r1=.1, repop.r2=.1, repop.s1=.1, repop.s2=.1, repop.s3=.1,
#                 short_dur=2, long_dur=10)

