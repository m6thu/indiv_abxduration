#########################################################################
#######Effect of antibiotic duration on hospitalised patients############
#########################################################################

rm(list=ls()) # Clean working environment
setwd("/Users/moyin/Desktop") #set working directory

###########################1. Input parameters##########################
#fixed parameters 
n.bed<-20                             # n.bed= number of beds in the ward
n.days<- 100                          # n.days= number of days we want to observe
mean.max.los<-10                      # mean.max.los= mean of max length of stay (normal distribution)
minDur<- 1                            # minDur=the min duration of antibiotic

#variable parameters 
###epidemiological 
p<-0.92                               # p=probability of receiving antibiotic
prob_StartBact<-c(0.4,0.58)           # prob_StartBact= vector of probability of carrying sensitive organism, resistant organism

###biological 
pi_s <- 0.003                         # pi_s= probability of S transmitting to N 
pi_r <- 0.003                         # pi_r= probability of R transmitting to N 
bif<- 1                               # bacterial interference factor 
pi_sr <- pi_r * bif                   # pi_sr= probability of R transmitting to S (a proportion of pi_r if being colonised with S protects colonisation by R)
mu_s <- 0                             # mu_s= rate of clearance of S to become N
mu_r <- 0                             # mu_r= rate of clearance of R to become S 

#########2. Generate length of stay and antibiotic duration table (1 for short duration and 1 for long duration) #######
#simulate inpatients with various lengths of stay
#allocate various duration of antibiotics for each patient

abx.table<- function (n.bed, n.days, mean.max.los, p, minDur, maxDur) {
  
  # generate a table of number of days we want to observe (rows) against number of beds in the ward (columns), filled in with patient id numbers
  patient.matrix<- matrix(NA, n.days, n.bed)  #make up a matrix of number of days we want to observe (rows) against number of beds in the ward (columns)
  n.patient=50*n.days                         #generate patient id numbers, arbitary number of 50 to ensure there are enough total number of patients generated to fill table
  n.patient.per.bed=c(1:n.patient)            #vectorise the patient id to be used for filling in the patient.matrix
  vector.los<-rep(0,n.days)                   #make up a fake vector for 0th column of matrix so that the next column starts with (last patient number of the previous column+i)
  
  max.los<- function() {                      # los is a normal function with mean of mean.max.los
    v<-round(rnorm(1, mean=mean.max.los, sd=2))
    if (v< 0) {
      v<--1*v
    }
    return(v)
  }
  
  for (j in 1:n.bed) {                                     #for each column (representing different beds in the ward) in patient.matrix
    los=c()
    final= vector.los[n.days]                              #last patient number in previous column
    
    for (i in 1:n.patient) {                               #for each row (representing number of days we want to observe)
      los <- rep(n.patient.per.bed[i+final], max.los())    #repeat the (last patient number of the previous column+i) random number of times
      vector.los <- c(vector.los, los)                     #combine 0th column (los0) and subsequent columns (los) into one vector
    }
    
    vector.los<- vector.los[-(1:n.days)]                   #remove 0th column from patient.matrix
    
    if (length(vector.los) > n.days) {                     #ensure each column is of the same length as number of days we are observing (number of rows)
      vector.los<-vector.los[1:n.days]
    }
    
    patient.matrix[,j] <- vector.los                       #fill the columns (different beds in the ward) with length of stay vector generated
  }
  
  #frequency summary of patient.matrix - patient id against number of days of stay for each patient
  los_duration<-as.vector(patient.matrix)
  
  #generate antibiotic use table
  
  # matrix_DuraDay is the matrix with number of days of antibiotics for every patient in patient.matrix
  matrix_DuraDay <- matrix(NA, nrow=n.days, ncol=n.bed)
  #returns first day of stay for patientnumber in bednumber, or NA if patient not there
  get1stdayofstay<-function(patientnumber,bednumber, bedoccmat){
    match(patientnumber, bedoccmat[,bednumber])
  }
  #returns los of patient patientnumber in bed bedumber, in bed occupancy matrix, 0 if patient not found
  getlos<-function(patientnumber,bednumber, bedoccmat){
    sum(bedoccmat[,bednumber]==patientnumber)
  }
  #number of days of antibiotic is randomly drawn from a uniform dist with min=minDur, max=maxDur
  for (i in 1:max(patient.matrix)){
    for (j in 1:n.bed){
      matrix_DuraDay[get1stdayofstay(i,j,patient.matrix), j]<-floor(runif(min=minDur,max=maxDur,1))
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
      matrix_AtbTrt[get1stdayofstay(i,j,patient.matrix), j] <-  if (rand <= p) {
        1
      }
      else{0}
    }
  }
  # Complete the matrix of cummulative length of stay for treated patients
  for (i in 2:n.days){
    for (j in 1:n.bed){
      matrix_AtbTrt[i, j] <-  if(is.na(matrix_AtbTrt[i,j]) & (matrix_AtbTrt[i-1,j]!=0)){
        matrix_AtbTrt[i-1,j] + 1
      }
      else if (!is.na(matrix_AtbTrt[i,j])){
        matrix_AtbTrt[i, j]
      }
      else{0}
    }
  }
  
  # Output matrix matrix_AtbTrt2 containing binary variable (treated vs not treated) for any bed on any particular day
  matrix_AtbTrt2 <- matrix(NA, nrow=n.days, ncol=n.bed)
  for (i in 1:n.days){
    for (j in 1:n.bed){
      matrix_AtbTrt2[i, j] <-  if(matrix_DuraDay[i,j]>=matrix_AtbTrt[i,j] & matrix_AtbTrt[i,j]!=0){
        1
      }
      else{0}
    }
  }
  
  return(list(patient.matrix, los_duration, matrix_AtbTrt2))
}

maxDur <-7                            # maxDur=the max duration of antibiotic
abx.short<-abx.table(n.bed, n.days, mean.max.los, p, minDur, maxDur)
maxDur <-14                            # maxDur=the max duration of antibiotic
abx.long<-abx.table(n.bed, n.days, mean.max.los, p, minDur, maxDur)

array_LOS<- function(los_duration) {
  los.dur<-as.vector(table(los_duration))
  array_LOS<-array(dim=c(2,length(los.dur)))
  array_LOS[1,]<-c(1:length(los.dur))
  array_LOS[2,]<-los.dur
  
  return(array_LOS)
}

array_LOS_short<-array_LOS(los_duration=abx.short[2])
array_LOS_long<- array_LOS(los_duration=abx.long[2])

#################3. Generate baseline carriage status ##################

gen_StartBact <- function(patient.matrix, prob_StartBact){
  #define probabilities of importing Sensitive(S) or Resistant(R) bacteria, or nothing (N)
  prob_start_S <- prob_StartBact[1]
  prob_start_R <- prob_StartBact[2]
  prob_start_N <- 1-prob_start_S-prob_start_R
  
  #Generating a vector of random status with runif (change for other distribution)
  number_of_patients<- max(as.data.frame(patient.matrix))
  Patient_StartBact <- runif(number_of_patients,0,1)
  Patient_StartBact[Patient_StartBact>prob_start_S+prob_start_R] <- 'N'
  Patient_StartBact[Patient_StartBact<=prob_start_S+prob_start_R&Patient_StartBact>prob_start_S] <- 'R'
  Patient_StartBact[Patient_StartBact<=prob_start_S] <- 'S'
  
  #Creating array for carriage status
  array_StartBact<-matrix(NA, n.days, n.bed)
  
  #returns first day of stay for patientnumber in bednumber, or NA if patient not there
  get1stdayofstay<-function(patientnumber,bednumber, bedoccmat){
    match(patientnumber, bedoccmat[,bednumber])
  }
  
  for (i in 1:number_of_patients){
    for (j in 1:n.bed){
      array_StartBact[get1stdayofstay(i,j, as.data.frame(patient.matrix)),j]<-Patient_StartBact[i]
    }
  }
  
  for (i in 1:n.days){
    for (j in 1:n.bed){
      array_StartBact[i, j] <-  if(is.na(array_StartBact[i,j])){
        0
      }
      else{array_StartBact[i, j]}
    }
  }
  
  return(array_StartBact)
}

array_StartBact_short<-gen_StartBact(patient.matrix=abx.short[1], prob_StartBact)
array_StartBact_long<-gen_StartBact(patient.matrix=abx.long[1], prob_StartBact)

####################4. Update values for every day  #####################
nextDay <- function(bed_table, array_LOS, treat_table, colo_table, pi_sr, mu_s, mu_r, pi_s, pi_r){
  
  # For each day (first day should be filled)
  for(i in 2:nrow(bed_table)){
    # For each bed
    for(j in 1:ncol(bed_table)){
      
      # case R
      #print(paste("i:", i, "j:", j))
      #print(colo_table[i-1, j])
      if(colo_table[i, j] == 0){
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
          if(treat_table[i-1, j] > 0){
            colo_table[i, j] <- "N"
          }else if(roll_clear < mu_s){
            # roll for clearance
            colo_table[i, j] <- "N"
          }else if(roll_select > pi_sr){
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
          prob_s <- 1-((1-.pi_s)^s_num)
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
      } # if 0
    }  # for j
  } # for i
  
  return(colo_table)
}

#output
colo_table_filled_short <- nextDay(bed_table= abx.short[[1]], array_LOS=array_LOS_short, treat_table=abx.short[[3]], colo_table=array_StartBact_short, pi_sr=pi_sr, mu_s=mu_s, mu_r=mu_r, pi_s=pi_s, pi_r=pi_r)
colo_table_filled_long <- nextDay(bed_table= abx.long[[1]], array_LOS=array_LOS_long, treat_table=abx.long[[3]], colo_table=array_StartBact_long, pi_sr=pi_sr, mu_s=mu_s, mu_r=mu_r, pi_s=pi_s, pi_r=pi_r)

####################6. Visualisation #####################

#present in one space
par(mfrow=c(4,2))

#plots
#antibiotics pixel
abx.mat.short<- as.matrix(abx.short[[3]])
abx.mat.long<- as.matrix(abx.long[[3]])
cols.a<-c('0'='grey', '1'='orange')
(abx_img.short<-image(1:nrow(abx.mat.short),1:ncol(abx.mat.short),abx.mat.short,col=cols.a, xlab='Time', ylab='Bed No.', main="Short duration of antibiotics"))
(abx_img.long<-image(1:nrow(abx.mat.long),1:ncol(abx.mat.long),abx.mat.long,col=cols.a, xlab='Time', ylab='Bed No.', main="Long duration of antibiotics"))

#carriage pixel
df.short<-as.data.frame(colo_table_filled_short)
df.long<-as.data.frame(colo_table_filled_long)
b.short<-data.matrix(df.short)
b.long<-data.matrix(df.long)
cols<-c('N'='grey', 'R'='red','S'='chartreuse3')
(col_img.short<-image(1:nrow(b.short),1:ncol(b.short),b.short,col=cols, xlab='Time', ylab='Bed No.'))
(col_img.long<-image(1:nrow(b.long),1:ncol(b.long),b.long,col=cols, xlab='Time', ylab='Bed No.'))

#increase overtime
df.short$totalR<-rep(0, nrow(df.short))
for (i in 1:nrow(df.short)) {
  for (j in 1:ncol(df.short)) {
    if (df.short[i,j]=="R") {
      df.short$totalR[i]<- df.short$totalR[i]+1
    } else {df.short$totalR[i]<- df.short$totalR[i]}
  }
}

df.long$totalR<-rep(0, nrow(df.long))
for (i in 1:nrow(df.long)) {
  for (j in 1:ncol(df.long)) {
    if (df.long[i,j]=="R") {
      df.long$totalR[i]<- df.long$totalR[i]+1
    } else {df.long$totalR[i]<- df.long$totalR[i]}
  }
}

x<-c(1:n.days)
y.short<-df.short$totalR/n.bed
y.long<-df.long$totalR/n.bed
r.plot.short<-plot(x,y.short, type="l", ylim=c(0,1), xlab="Time", ylab="Proportion of patients carrying resistant organisms")
r.plot.long<-plot(x, y.long, type="l", ylim=c(0,1), xlab="Time", ylab="Proportion of patients carrying resistant organisms")

#total acquisitions 
df.short$newR<-rep(0, nrow(df.short))
df.long$newR<-rep(0, nrow(df.long))

for(i in 2:nrow(df.short)){
    for(j in 1:ncol(df.short)){
      if(df.short[i, j] == "R"){
        if(df.short[i-1, j] == "S"|df.short[i-1, j] == "N"){
          df.short$newR[i]<-df.short$newR[i]+1 
      } else {df.short$newR[i]<-df.short$newR[i]}
    }
  }
}
y.short<- df.short$newR

for(i in 2:nrow(df.long)){
  for(j in 1:ncol(df.long)){
    if(df.long[i, j] == "R"){
      if(df.long[i-1, j] == "S"|df.long[i-1, j] == "N"){
        df.long$newR[i]<-df.long$newR[i]+1 
      } else {df.long$newR[i]<-df.long$newR[i]}
    }
  }
}
y.long<- df.long$newR

acq.short<-plot(x,y.short, type='l', xlab="Time", ylab="Number of patients acquiring resistant organisms")
acq.long<-plot(x,y.long, type='l', xlab="Time", ylab="Number of patients acquiring resistant organisms")
