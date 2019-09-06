#######Modelling Day Project######
#######Parameter exploration######
##################################

# SAMPLE PARAMETER SPACE 
# load libraries 
library(pse) #load pse package for Latin Hypercube
library (sensitivity) #load sensitivity package for sensitivity analysis 

#binary model functions 
######## Functions ########
positive_norm_sample <- function(mean, sd){                      
    v<-round(rnorm(1, mean=mean, sd=sd))
    if (v < 0) {
        v<-v*-1
    }
    return(v)
}
#########2. Generate length of stay and antibiotic duration table
#(1 for short duration and 1 for long duration) 
#simulate inpatients with various lengths of stay
#allocate various duration of antibiotics for each patient

abx.table<- function (n.bed, mean.max.los, p.s, p.r, meanDur) {
    
    # generate a table of number of days we want to observe (rows) -
    # against number of beds in the ward (columns), filled in with patient id numbers
    
    n.days=200
    
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
    
    # los is a normal function with mean of mean.max.los
    # positive_norm_sample
    
    for (j in 1:n.bed) {
        #for each column (representing different beds in the ward) in patient.matrix
        
        los <- c()
        final <- vector.los[n.days]
        #last patient number in previous column
        
        for (i in 1:n.days) {
            #for each row (representing number of days we want to observe)
            
            los <- rep(patient.id[i+final], positive_norm_sample(mean.max.los, sd=2))
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
            if (rand < p.r) {
                matrix_AtbTrt.r[get1stdayofstay(i,j,patient.matrix), j] <-  2
            } else {
                matrix_AtbTrt.r[get1stdayofstay(i,j,patient.matrix), j] <-  0
            }
        }
    }
    # Initial treatment value derived from probability, p.r for antibiotic.r 
    
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
    # howmanydaysofabt(m, 2, 1) should return 1
    # howmanydaysofabt(m, 2, 2) should return 1 < still problematic
    # Check cross over between patients
    # howmanydaysofabt(m, 4, 1) should return 2
    # howmanydaysofabt(m, 4, 3) should return 3
    
    
    for (i in 2:nrow(patient.matrix)){
        for (j in 1:n.bed){
            #print(paste("i", i, "j", j))
            #print(matrix_AtbTrt.r[i, j])
            if(is.na(matrix_AtbTrt.r[i, j])){
                rand <- runif(1,0,1)
                # case of no antibiotics for resistant organisms in the day before
                if (matrix_AtbTrt.r[i-1, j] == 0) {
                    if (rand < p.r) {
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

gen_StartBact <- function(los, n.bed, prob_StartBact_R){
    
    n.days=200
    
    prob_start_S <-  0.75*(1- prob_StartBact_R)
    prob_start_s <-  0.25*(1- prob_StartBact_R)
    prob_start_Sr <-  0.8* prob_StartBact_R
    prob_start_sr <-  0.1* prob_StartBact_R
    prob_start_sR <-  0.1* prob_StartBact_R
    prob_StartBact_bi<- c(prob_start_S,prob_start_Sr,prob_start_sR,prob_start_sr) 
    
    stopifnot(sum(prob_StartBact_bi) < 1) # Assert all probabilities combined are less than 1
    
    #define probabilities of importing S, Sr, sR, sr, s
    # prob_start_S <-  prob_StartBact[1]
    # prob_start_Sr <- prob_StartBact[2]
    # prob_start_sR <- prob_StartBact[3]
    # prob_start_sr <- prob_StartBact[4]
    # prob_start_s <- 1 - sum(prob_StartBact)
    
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
    
    return(array_StartBact)
}


####################4. Update values for every day  
nextDay <- function(bed_table, array_LOS, treat_table, colo_table, pi_r1, pi_r2, mu1, mu2, repop.r1, repop.r2, repop.s1, repop.s2,repop.s3, abx.r,abx.s){
    
    # For each day (first day should be filled)
    for(i in 2:nrow(bed_table)){
        # For each bed
        for(j in 1:ncol(bed_table)){
            #case S
            #print(paste("i:", i, "j:", j))
            #print(colo_table[i-1, j])
            if(is.na(colo_table[i, j])){
                if(colo_table[i-1, j] == "S"){
                    #print("----case S")
                    # check antibiotic
                    roll_clear <- runif(1, 0, 1)
                    roll_transmit <- runif(1, 0, 1)
                    if (treat_table[i-1, j] == 1 & roll_clear < abx.s){
                        colo_table[i, j] <- "s"
                    } else if (treat_table[i-1, j] > 1 & roll_clear < abx.r){
                        colo_table[i, j] <- "s"
                    } else if (roll_transmit < pi_r1){ 
                        colo_table[i, j] <- "Sr"
                    }else {
                        colo_table[i, j] <- "S"
                    }
                    
                    # case s
                }else if(colo_table[i-1, j] == "s"){
                    #print("----case s")
                    # roll for transmission of r
                    roll_r <- runif(1, 0, 1)
                    r_num <- sum(colo_table[i-1,] == "sR") #only R can be transmitted 
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

#list parameters, the probability density functions from which the parameter values will be calculated, and what are the arguments to these functions
factors <- c(             #list parameters in an array 
    "n.bed",              #number of beds in the ward
    #"n.days",             #number of days of observation 
    "mean.max.los",       #mean of length of stay (normal distribution)
    "p.s",                #probability of being prescribed narrow spectrum antibiotic
    "p.r",                #probability of being prescribed broad spectrum antibiotic
    "prob_StartBact_R",   #probability of initial carriage of resistant organisms
    "pi_r1",              #probability of being transmitted r to S (S—> Sr)
    "pi_r2",              #probability of being transmitted r to s (s—>sr)
    "mu1",                #probability of being decolonised to S (Sr—> S) 
    "mu2",                #probability of being decolonised to S (sr—> s) 
    "abx.r",              #probability of clearing R to become r
    "abx.s",              #probability of clearing S to become s
    "repop.r1",           #probability of transmission of r to S (s—> Sr) 
    "repop.r2",           #probability of regrowth of s (sR—> sr)
    "repop.s1",           #probability of regrowth of S  (s—>S)
    "repop.s2",           #probability of regrowth of S  (sr—>Sr)
    "repop.s3",           #probability of transmission of r to S (s—> Sr) 
    #"bif",               #bacterial interference factor (pi_r2 = pi_r1 * bif )
    #"bif1",              #bacterial interference factor 1 (repop.r3=repop.r2*bif1)
    "short_dur",          #mean short duration of antibiotics (normal distribution) 
    "long_dur"            #mean long duration of antibiotics (normal distribution) 
    )  
q <- c(                   #distributions of parameters 
    "qunif",              # "n.bed", number of beds in the ward
    #"qunif",              # "n.days" number of days of observation 
    "qunif",              # "mean.max.los", mean of length of stay (normal distribution)
    "qunif",              # "p.s", probability of being prescribed narrow spectrum antibiotic
    "qunif",              # "p.r", probability of being prescribed broad spectrum antibiotic
    "qunif",              # "prob_StartBact_R",probability of initial carriage of resistant organisms
    "qunif",              # "pi_r1" probability of being transmitted r to S (S—> Sr)
    "qunif",              # "pi_r2",probability of being transmitted r to s (s—>sr)
    "qunif",              # "mu1" probability of being decolonised to S (Sr—> S) 
    "qunif",              # "mu2",probability of being decolonised to S (sr—> s) 
    "qunif",              # "abx.r", probability of clearing R to become r
    "qunif",              # "abx.s", probability of clearing S to become s
    "qunif",              # "repop.r1", probability of transmission of r to S (s—> Sr) 
    "qunif",              # "repop.r2", probability of regrowth of s (sR—> sr)
    "qunif",              # "repop.s1", probability of regrowth of S  (s—>S)
    "qunif",              # "repop.s2", probability of regrowth of S  (sr—>Sr)
    "qunif",              # "repop.s3",probability of transmission of r to S (s—> Sr) 
    "qunif",              #"short_dur", mean short duration of antibiotics (normal distribution) 
    "qunif"               #"long_dur", mean long duration of antibiotics (normal distribution)  
    ) 

q.arg <- list(            #set limits of parameters 
    list(min=10, max=20),                # "n.bed", number of beds in the ward
    #list(min=10, max=20),                # "n.days" number of days of observation 
    list(min=5, max=10),                 # "mean.max.los", mean of length of stay (normal distribution)
    list(min=0.1, max=0.9),              # "p.s", probability of being prescribed narrow spectrum antibiotic
    list(min=0.1, max=0.9),              # "p.r", probability of being prescribed broad spectrum antibiotic
    list(min=0.1, max=0.7),              # "prob_StartBact_R",probability of initial carriage of resistant organisms
    list(min=0.1, max=0.9),              # "pi_r1" probability of being transmitted r to S (S—> Sr)
    list(min=0.1, max=0.9),              # "pi_r2",probability of being transmitted r to s (s—>sr)
    list(min=0.1, max=0.9),              # "mu1" probability of being decolonised to S (Sr—> S) 
    list(min=0.1, max=0.9),              # "mu2",probability of being decolonised to S (sr—> s) 
    list(min=0.1, max=0.9),              # "abx.r", probability of clearing R to become r
    list(min=0.1, max=0.9),              # "abx.s", probability of clearing S to become s
    list(min=0.1, max=0.9),              # "repop.r1", probability of transmission of r to S (s—> Sr) 
    list(min=0.1, max=0.9),              # "repop.r2", probability of regrowth of s (sR—> sr)
    list(min=0.1, max=0.9),              # "repop.s1", probability of regrowth of S  (s—>S)
    list(min=0.1, max=0.9),              # "repop.s2", probability of regrowth of S  (sr—>Sr)
    list(min=0.1, max=0.9),              # "repop.s3",probability of transmission of r to S (s—> Sr) 
    list(min=4, max=10),                 #"short_dur", mean short duration of antibiotics (normal distribution) 
    list(min=14, max=20)                 #"long_dur", mean long duration of antibiotics (normal distribution)  
    )

#model 

#one run of the model 
diff_prevalence <- function(n.bed, mean.max.los, p.s, p.r,
                            prob_StartBact_R, pi_r1, pi_r2, mu1, mu2, abx.r, abx.s,
                            repop.r1, repop.r2, repop.s1, repop.s2, repop.s3,
                            short_dur, long_dur){
    n.days=200
    iterations=10
    iter_totalsR <- matrix(NA, nrow = n.days, ncol = iterations)
    
    for(iter in 1:iterations){
        
        #print(paste("iter:", iter, "y:", y_count, '-', y, "x", x_count, '-', x))
        #Generate length of stay and antibiotic duration table
        abx_iter <- abx.table(n.bed=n.bed, mean.max.los=mean.max.los, p.s=p.s, p.r=p.r, meanDur=short_dur)
        #Generate baseline carriage status
        array_LOS_iter <- array_LOS_func(los_duration=abx_iter[[2]])
        #Update values for every day
        array_StartBact_iter <- gen_StartBact(los=array_LOS_iter, n.bed=n.bed, prob_StartBact_R=prob_StartBact_R)
        #output
        colo_table_filled_iter <- nextDay(bed_table= abx_iter[[1]], array_LOS=array_LOS_iter, 
                                          treat_table=abx_iter[[3]], colo_table=array_StartBact_iter, 
                                          pi_r1=pi_r1, pi_r2= pi_r2, mu1=mu1, mu2=mu2, 
                                          abx.r=abx.r,abx.s=abx.s,
                                          repop.r1 = repop.r1, repop.r2 = repop.r2,
                                          repop.s1 = repop.s1, repop.s2 = repop.s2,repop.s3 = repop.s3)
        #Summary
        df <- data.frame(colo_table_filled_iter)
        iter_totalsR[, iter] <- rowSums(df == "sR")
        #print("end iteration loop")
    }
    totalsR_short <- mean(rowSums(iter_totalsR)/iterations/n.bed)
    
    iter_totalsR <- matrix(NA, nrow = n.days, ncol = iterations)
    for(iter in 1:iterations){
        
        #print(paste("iter:", iter, "y:", y_count, '-', y, "x", x_count, '-', x))
        #Generate length of stay and antibiotic duration table
        abx_iter <- abx.table(n.bed=n.bed, mean.max.los=mean.max.los, p.s=p.s, p.r=p.r, meanDur=long_dur)
        #Generate baseline carriage status
        array_LOS_iter <- array_LOS_func(los_duration=abx_iter[[2]])
        #Update values for every day
        array_StartBact_iter <- gen_StartBact(los=array_LOS_iter, n.bed=n.bed,prob_StartBact_R=prob_StartBact_R)
        #output
        colo_table_filled_iter <- nextDay(bed_table= abx_iter[[1]], array_LOS=array_LOS_iter, 
                                          treat_table=abx_iter[[3]], colo_table=array_StartBact_iter, 
                                          pi_r1=pi_r1, pi_r2= pi_r2, mu1=mu1, mu2=mu2, 
                                          abx.r=abx.r,abx.s=abx.s,
                                          repop.r1 = repop.r1, repop.r2 = repop.r2,
                                          repop.s1 = repop.s1, repop.s2 = repop.s2, repop.s3 = repop.s3)
        #Summary
        df <- data.frame(colo_table_filled_iter)
        iter_totalsR[,iter] <- rowSums(df == "sR")
        #print("end iteration loop")
    }
    totalsR_long <- mean(rowSums(iter_totalsR)/iterations/n.bed)
    
    print(paste("totalsR_long", totalsR_long, "totalsR_short", totalsR_short))
    
    return(totalsR_long - totalsR_short)
}

modelRun.binary <- function (data.df) { #data.df is a dataframe of the parameter values in columns 
    return(mapply(diff_prevalence, 
                  data.df[,1], data.df[,2], data.df[,3], 
                  data.df[,4], data.df[,5], data.df[,6], 
                  data.df[,7], data.df[,8], data.df[,9], 
                  data.df[,10], data.df[,11], data.df[,12], 
                  data.df[,13], data.df[,14], data.df[,15], 
                  data.df[,16], data.df[,17], data.df[,18]
                  ))
    }

# Use the LHD function to generate a hypercube 
LHS.binary<- LHS(modelRun.binary, factors, 20, q, q.arg, nboot=20)
results.binary<-get.results(LHS.binary)

#Plot findings: 
#1. empirical cumulative distribution function used to illustrate the distribution of the model results
plotecdf(LHS.binary) #outcome has a high probability in the steepest parts of the graph 

#2. scatterplot of the result as a function of each parameter: distribution of values returned by the model 
# in the parameter space sampled by the hypercube and how sensible are these model responses to the variation of each parameter.
plotscatter(LHS.binary)

#3. partial correlation coefficient measures how strong are the inear associations between the result 
# and each input parameter, after removing the linear effect of the other parameters. (CI generated by boostrapping)
plotprcc(LHS.binary)

pic(LHS.binary, nboot=40) #pic (partial inclination coefficient) is the sensitivity" of the model response in respect to each parameter
#represent the beta terms in y = alpha + beta*x regressions, after removing the linear effect of the other parameters

# 4. Cobweb
outcome.df<-as.data.frame(cbind(LHS.binary$data,results.binary)) #dummy matrix with parameter values in columns and outcome in last column
names(outcome.df)<- c(factors,'outcome') #name the columns of the dummy matrix 
for (i in 1:nrow(outcome.df)) {       #label the rows of parameter values that produced top 5% of the outcomes
    if (outcome.df$outcome[i]<quantile(outcome.df$outcome,probs = 0.95)) { 
        outcome.df$top5[i] <-0 } else {
            outcome.df$top5[i] <-1
        }
}
library(MASS) #load MASS package
colors<- c("#E69F00", "#009E73") #choose 2 colors - 1 for parameters that produced top 5% of outcomes and one for the rest
outcome.df$top5<- as.factor(outcome.df$top5)
parcoord(outcome.df[,c(1:length(factors))] , col= colors[outcome.df$top5],var.label=T)

#5. Check agreement between runs to decide if our sample size for adequate 
# Symmetric Blest Measure of Agreement (SBMA) between the PRCC coeffients of two runs with different sample sizes.
check.LHS.binary <- LHS(modelRun.binary, factors.binary, 250, q.binary, q.arg.binary)
(mySbma <- sbma(LHS.binary, check.LHS.binary))
# value of -1 indicates complete disagreement between the runs 
# value of 1 indicated complete agreement  (>0.7 acceptable)
#caveat: if none of the model parameters is monotonically correlated with the output, the agreement between runs may stay as low as 0.2 even for very large hypercubes.

#################################################################################
###DUMMY CODE FOR PARAMETER EXPLORATION###
fun.trial<- function(x1, x2, x3) {x1+x2^2+ x3^3} #dummy model 
factors.trial<-c('x1', 'x2', 'x3') #dummy variables 
q.trial<- c("qunif","qunif","qunif") #distribution of the variables 
q.arg.trial <- list(        #set limits to the variables   
    list(min=0.1, max=0.9),     
    list(min=0.1, max=0.9),       
    list(min=0.1, max=0.9)) 
modelRun.trial <- function (data.df) { #data.df is a dataframe of the parameter values in columns 
    return(mapply(fun.trial, data.df[,1], data.df[,2], data.df[,3]))
}

LHS.trial<- LHS(modelRun.trial, factors.trial, 50, q.trial, q.arg.trial, nboot=20) #50 parameter combinations to be generated, nboot= number of correlation coefficients 
results.trial<-get.results(LHS.trial)

#Plot findings: 
#1. empirical cumulative distribution function used to illustrate the distribution of the model results
plotecdf(LHS.trial) #outcome has a high probability in the steepest parts of the graph 

#2. scatterplot of the result as a function of each parameter: distribution of values returned by the model 
# in the parameter space sampled by the hypercube and how sensible are these model responses to the variation of each parameter.
plotscatter(LHS.trial)

#3. partial correlation coefficient measures how strong are the inear associations between the result 
# and each input parameter, after removing the linear effect of the other parameters. (CI generated by boostrapping)
plotprcc(LHS.trial)

pic(LHS.trial, nboot=40) #pic (partial inclination coefficient) is the sensitivity" of the model response in respect to each parameter
#represent the beta terms in y = alpha + beta*x regressions, after removing the linear effect of the other parameters

# 4. Cobweb
outcome.df<-as.data.frame(cbind(LHS.trial$data,results.trial)) #dummy matrix with parameter values in columns and outcome in last column
names(outcome.df)<- c('X1','X2','X3','Y') #name the columns of the dummy matrix 
for (i in 1:nrow(outcome.df)) {       #label the rows of parameter values that produced top 5% of the outcomes
    if (outcome.df$Y[i]<quantile(outcome.df$Y,probs = 0.95)) { 
        outcome.df$top5[i] <-0 } else {
            outcome.df$top5[i] <-1
        }
}
library(MASS) #load MASS package
colors<- c("#E69F00", "#009E73") #choose 2 colors - 1 for parameters that produced top 5% of outcomes and one for the rest
outcome.df$top5<- as.factor(outcome.df$top5)
parcoord(outcome.df[,c(1:3)] , col= colors[outcome.df$top5],var.label=T)

#5. Check agreement between runs to decide if our sample size for adequate 
# Symmetric Blest Measure of Agreement (SBMA) between the PRCC coeffients of two runs with different sample sizes.
check.LHS.trial <- LHS(modelRun.trial, factors.trial, 250, q.trial, q.arg.trial)
(mySbma <- sbma(LHS.trial, check.LHS.trial))
# value of -1 indicates complete disagreement between the runs 
# value of 1 indicated complete agreement  (>0.7 acceptable)
#caveat: if none of the model parameters is monotonically correlated with the output, the agreement between runs may stay as low as 0.2 even for very large hypercubes.
