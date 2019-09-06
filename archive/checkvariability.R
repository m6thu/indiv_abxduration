# Adapted from MDP_individualmodel_110918.R

#########################################################################
#######Effect of antibiotic duration on hospitalised patients############
#########################################################################

rm(list=ls()) # Clean working environment
setwd("~/Desktop") #set working directory

###########################1. Input parameters##########################
#fixed parameters 
n.bed<-20                             # n.bed= number of beds in the ward
n.days<- 100                          # n.days= number of days we want to observe
mean.max.los<-20                      # mean.max.los= mean of max length of stay (normal distribution)

#variable parameters 
###epidemiological 
p.s<-0                              # p=probability of receiving antibiotic for sensitive organisms
p.r<-0                              # p= daily probability of contracting HAI and receiving antibiotic for resistant organisms 
prob_StartBact<-c(0.5,0.2, 0.1, 0.05) # prob_StartBact= vector of probability of carrying c(S,Sr, sR, sr)
#                                     # possible states: S- carry sensitive organism only 
#                                                        Sr- carry largely sensitive organisms and small population of resistant organisms
#                                                        sR- carry largely resistant organisms and small population of sensitive organisms
#                                                        sr- carry small population of sensitive organisms and resistant organisms
#                                                        s - carry small population of sensitive organisms 

###biological 
pi_r1 <- 0.003                        # pi_r1= probability of R transmitting to S to become Sr
bif<- 2                               # bacterial interference factor 
pi_r2 <- pi_r1 * bif                  # pi_r2= probability of R transmitting to s to become sr 
#                                       (pi_r1 < pi_r2 if being colonised with S protects colonisation by R)

mu1 <- 0                              # mu1= probability of clearance of Sr to become S
mu2 <- 0                              # mu2= probability of clearance of sr to become s 

abx.s<-0.2                            # probability of clearing S to become s under antibiotic treatment (daily)
abx.r<-0.3                            # probability of clearing R to become r under antibiotic treatment (daily)

repop.s1<- 0                          # probability of repopulation of s to become S 
repop.s2<- 0                          # probability of repopulation of sr to become SR 
repop.s3<- 0                       # probability of repopulation of sR to become sr
repop.r1<- 0                          # probability of repopulation of Sr to become sR 
repop.r2<- 0                          # probability of repopulation of sr to become sR 

positive_norm_sample <- function(mean, sd){                      
    v<-round(rnorm(1, mean=mean, sd=sd))
    if (v < 0) {
        v<-v*-1
    }
    return(v)
}

###################### Define Functions #############################

#2. Generate length of stay and antibiotic duration table (define function)
#(1 for short duration and 1 for long duration) 
#simulate inpatients with various lengths of stay
#allocate various duration of antibiotics for each patient

abx.table<- function (n.bed, n.days, mean.max.los, p.s, p.r, meanDur) {
    
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

array_LOS<- function(los_duration) {
    los.dur<-as.vector(table(los_duration))
    array_LOS<-array(dim=c(2,length(los.dur)))
    array_LOS[1,]<-c(1:length(los.dur))
    array_LOS[2,]<-los.dur
    
    return(array_LOS)
}

#3. Generate baseline carriage status (define function)

gen_StartBact <- function(los, prob_StartBact){
    
    stopifnot(sum(prob_StartBact) < 1) # Assert all probabilities combined are less than 1
    
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
    Patient_StartBact[Patient_unif > sum(prob_StartBact)] <- 's'
    Patient_StartBact[(Patient_unif > sum(prob_StartBact[1:3])) & (Patient_unif <= sum(prob_StartBact))] <- 'sr'
    Patient_StartBact[(Patient_unif > sum(prob_StartBact[1:2])) & (Patient_unif <= sum(prob_StartBact[1:3]))] <- 'sR'
    Patient_StartBact[(Patient_unif > sum(prob_StartBact[1])) & (Patient_unif <= sum(prob_StartBact[1:2]))] <- 'Sr'
    Patient_StartBact[Patient_unif <= sum(prob_StartBact[1])] <- 'S'
    
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

# 4. Update values for every day (define function)
nextDay <- function(bed_table, array_LOS, treat_table, colo_table, pi_r1, pi_r2, mu1, mu2, repop.r1, repop.r2, repop.s1, repop.s2){
    
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
                    }else if(treat_table[i-1, j] == 0 &roll_repop < repop.s2){
                        colo_table[i, j] <- "Sr"
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

################################## Multiple runs comparing between both ####################################
iterations <- 100

abx.short <- list()
abx.long <- list()
array_LOS_short <- list()
array_LOS_long <- list()
array_StartBact_short <- list()
array_StartBact_long <- list()
colo_table_filled_short <- list()
colo_table_filled_long <- list()

short_totalR <- matrix(NA, nrow = n.days, ncol = iterations)
long_totalR <- matrix(NA, nrow = n.days, ncol = iterations)

mean_short <- 4 #antibiotic short duration
mean_long <- 14 #antibiotic short duration

for (i in 1:iterations){
    
    print(i)
    #Generate length of stay and antibiotic duration table
    abx.short[[i]] <- abx.table(n.bed, n.days, mean.max.los, p.s, p.r, meanDur=mean_short)
    abx.long[[i]] <- abx.table(n.bed, n.days, mean.max.los, p.s, p.r, meanDur=mean_long)
    
    #Generate baseline carriage status
    array_LOS_short[[i]] <- array_LOS(los_duration=abx.short[[i]][2])
    array_LOS_long[[i]] <- array_LOS(los_duration=abx.long[[i]][2])
    
    #Update values for every day
    array_StartBact_short[[i]] <- gen_StartBact(los=array_LOS_short[[i]], prob_StartBact)
    array_StartBact_long[[i]] <- gen_StartBact(los=array_LOS_long[[i]], prob_StartBact)
    
    #output
    colo_table_filled_short[[i]] <- nextDay(bed_table= abx.short[[i]][[1]], array_LOS=array_LOS_short[[i]], 
                                            treat_table=abx.short[[i]][[3]], colo_table=array_StartBact_short[[i]], 
                                            pi_r1=pi_r1, mu1=mu1, mu2=mu2, pi_r2=pi_r2, repop.r1 = repop.r1, 
                                            repop.r2 = repop.r2, repop.s1 = repop.s1, repop.s2 = repop.s2)
    colo_table_filled_long[[i]] <- nextDay(bed_table= abx.long[[i]][[1]], array_LOS=array_LOS_long[[i]], 
                                           treat_table=abx.long[[i]][[3]], colo_table=array_StartBact_long[[i]], 
                                           pi_r1=pi_r1, mu1=mu1, mu2=mu2, pi_r2=pi_r2, repop.r1 = repop.r1, 
                                           repop.r2 = repop.r2, repop.s1 = repop.s1, repop.s2 = repop.s2)

    #increase overtime
    df.short <- colo_table_filled_short[[i]]
    short_totalR[,i]<- rowSums(df.short == "sR")
    
    df.long <- colo_table_filled_long[[i]]
    long_totalR[,i] <- rowSums(df.long == "sR")

}

# Debug edge cases: 43, 71 from run 1431_28Sept2018.RData
# Debug cases: 36, 61, 160, 163, 212, 240, 258, 275, 300, 557, 565, 648, 658, 692, 735, 783, 813, 897 from run 1501_28Sept2018.RData
# A clean run: 1553_28Sept2018.RData, still uses lists in iteration collection instead of matrix that is easier to summarize from


# Check between 2 simple scenarios
par(mfrow=c(1,2))

x <- 1:n.days
y.short <- short_totalR[,1]/n.bed
r.plot.short <- plot(x, y.short, type="l", ylim=c(0,1), main=paste("Short mean abx:", mean_short, "days", "(", round(mean(rowSums(short_totalR)/iterations/n.bed), digits=2), ")"),
                     xlab="Time", ylab="% patients carrying MDRO", col=rgb(0,0,0,alpha=0.1))
for(i in 2:iterations){
    lines(x, short_totalR[,i]/n.bed, col=rgb(0,0,0,alpha=0.1))
}
lines(x, rowSums(short_totalR)/iterations/n.bed, col=2)
abline(h = mean(rowSums(short_totalR)/iterations/n.bed), col=2, lwd=2)

y.long <- long_totalR[,1]/n.bed
r.plot.long <- plot(x, y.long, type="l", ylim=c(0,1), main=paste("Long mean abx:", mean_long, "days", "(", round(mean(rowSums(long_totalR)/iterations/n.bed), digits=2), ")"),
                    xlab="Time", ylab="% patients carrying MDRO", col=rgb(0,0,0,alpha=0.1))
for(i in 2:iterations){
    lines(x, long_totalR[,i]/n.bed, col=rgb(0,0,0,alpha=0.1))
}
lines(x, rowSums(long_totalR)/iterations/n.bed, col=2)
abline(h = mean(rowSums(long_totalR)/iterations/n.bed), col=2, lwd=2)
