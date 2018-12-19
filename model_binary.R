require(msm)

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


# Rename to replace bottom (comment out) when its been tested to work
abx.table <- function(patient.matrix, los.array, p.s, p.r.day1, p.r.dayafter,
                          meanDur.s, meanDur.r, sdDur, timestep=1){

    # Check assumption that possibilities are, no abx, has s abx, or has r abx on first day
    stopifnot(p.s+p.r.day1 <= 1)

    #generate antibiotic use table
    #number of days of s antibiotic is randomly drawn from a truncated normal distribution
    abx_days.s <- round(rtnorm(ncol(los.array), mean=meanDur.s*timestep, sd=sdDur*timestep, lower=1))
    #number of days of r antibiotic is randomly drawn from a truncated normal distribution
    abx_days.r <- round(rtnorm(ncol(los.array), mean=meanDur.r*timestep, sd=sdDur*timestep, lower=1))
    #number of days of r antibiotic for days after drawn from tunced norm dist
    abx_r.after <- round(rtnorm(ncol(los.array), mean=meanDur.r*timestep, sd=sdDur*timestep, lower=1))
    r_idx <- 1 # R indices start at 1
    # Unit test - check distribution of abx distribution
    # hist(abx_days.s, breaks=20)
    # Unit test - compare cases that will enter padding if-else
    # abx_days > los.array[2, ]

    # abx.matrix should be same size as patient.matrix
    abx.matrix <- matrix(NA, nrow=nrow(patient.matrix), ncol=ncol(patient.matrix))
    idx_end <- 1
    # For each patient
    for(i in 1:ncol(los.array)){

        # maxiumum number of days for a particular patient from los vector
        max_days <- los.array[2, i]
        # number of abx days for that particular patient from generated number
        abx_s <- abx_days.s[i]
        abx_r <- abx_days.r[i]
        
        # Initial treatment value derived from probability for p.r and p.s
        rand <- runif(1, 0, 1)
        if (rand < p.s){ # check if p.s starting
            if(abx_s > max_days){
                # if the number of generated abx is longer than the los of that person
                # have them take abx everyday for their stay
                abx.matrix[idx_end:(idx_end+max_days-1)] <- rep(1, max_days)
            }else{
                # else take abx for abx days and pad to fit max_days
                abx.matrix[idx_end:(idx_end+max_days-1)] <- c(rep(1, abx_s), rep(0, max_days-abx_s))
            }
        }else if(rand < (p.s+p.r.day1)){ # case for p.r starting, rand >= p.s already checked in previous if case
            if(abx_r > max_days){
                # if the number of generated abx is longer than the los of that person
                # have them take abx everyday for their stay
                abx.matrix[idx_end:(idx_end+max_days-1)] <- rep(2, max_days)
            }else{
                # else take abx for abx days and pad to fit max_days
                abx.matrix[idx_end:(idx_end+max_days-1)] <- c(rep(2, abx_r), rep(0, max_days-abx_r))
            }
        }else{
            # no abx taken for that person
            abx.matrix[idx_end:(idx_end+max_days-1)] <- rep(0, max_days)
        }
        
        #Every day has a chance of starting abx.r by a fixed probability
        roll.r <- runif(max_days, 0, 1)
        start.r <- roll.r < p.r.dayafter
        #print(sum(start.r))
        where.r <- which(start.r == 1)
        #print(paste('where.r', where.r))
        # Create an array of p.r starting based on this
        # Merge starting vector and random p.r start vector (giving priority to p.r random start)
        # for each location that rolled chance that abx.r starts
        if(sum(start.r)){
            for(j in 1:sum(start.r)){
                #print(j)
                # replace that location for the length of abx_r.after drawn from norm distribution
                start_idx <- idx_end+where.r[j]-1
                #print(start_idx)
                end_idx <- start_idx+abx_r.after[r_idx]-1
                #print(end_idx)
                if(end_idx > (idx_end+max_days-1)){ # abx duration exceeds max_days for that person, replace the rest, escape loop
                    abx.matrix[start_idx:(idx_end+max_days-1)] <- rep(2, max_days-where.r[j]+1)
                    break
                }else{ # abx duration does not exceed max_days for that person
                    abx.matrix[start_idx:end_idx] <- rep(2, abx_r.after[r_idx])
                }
                # move to next abx_r.after drawn
                r_idx <- r_idx+1
            }
        }
        
        # move starting position to end of previous patient
        idx_end <- idx_end+max_days
    }

    return(abx.matrix)
}



abx.table.old <- function(patient.matrix, los.array, p.s, p.r.day1, p.r.dayafter,
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
nextDay.old <- function(patient.matrix, abx.matrix, colo.matrix, 
                    pi_r1, bif, mu1, mu2, repop.r1, repop.r2, 
                    repop.s1, repop.s2,depop.r, abx.r, abx.s, timestep=1){
    
    # adjust probabilities based on timestep
    pi_r1 <- pi_r1/timestep
    mu1 <- mu1/timestep
    mu2 <- mu2/timestep
    repop.r1 <- repop.r1/timestep
    repop.r2 <- repop.r2/timestep
    repop.s1 <- repop.s1/timestep
    repop.s2 <- repop.s2/timestep
    depop.r <- depop.r/timestep
    abx.r <- abx.r/timestep
    abx.s <- abx.s/timestep
    
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
                    }else if (abx.matrix[i-1, j] == 0 & roll_clear < depop.r){
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

nextDay <- function(patient.matrix, abx.matrix, colo.matrix, 
                        pi_r1, bif, mu1, mu2, repop.r1, repop.r2, 
                        repop.s1, repop.s2, depop.r, abx.r, abx.s, timestep=1){
    
    # adjust probabilities based on timestep
    pi_r1 <- pi_r1/timestep
    mu1 <- mu1/timestep
    mu2 <- mu2/timestep
    repop.r1 <- repop.r1/timestep
    repop.r2 <- repop.r2/timestep
    repop.s1 <- repop.s1/timestep
    repop.s2 <- repop.s2/timestep
    depop.r <- depop.r/timestep
    abx.r <- abx.r/timestep
    abx.s <- abx.s/timestep
    
    pi_r2 <- pi_r1 * bif                 # pi_r2= probability of R transmitting to s to become sr 
    
    # For each day (first day should be filled)
    for(i in 2:nrow(patient.matrix)){
        # Get the previous row subset of the entire matrix which represents previous day or timestep
        prev_step <- colo.matrix[i-1, ]
        # Get the indices which are already filled by a patient entering the ward for exclusion
        already_filled <- which(!is.na(colo.matrix[i, ]))
        # only R can be transmitted, get number of r for transmission probability
        #print(colo.matrix[i-1,])
        r_num <- sum(colo.matrix[i-1,] == "sR") 
        
        # Update S
        # Get the column indices which contain S in the previous day
        S <- which(prev_step == "S")
        # Remove column indices that already have a starting bacterial state filled in
        S <- S[!(S %in% already_filled)]
        # count if there are any S on the previous day
        s_num <- length(S)
        # if there is any S (number of S > 0) in the previous day
        if(s_num){
            # roll for transmission of R
            prob_r <- 1-((1-pi_r1)^r_num)
            # Roll a random number for each R on the previous day for clearance
            roll_clear <- runif(s_num, 0, 1)
            # All column indices which roll < probability of clearance AND there is antibiotic in the previous day
            clear_idx.s <- S[abx.matrix[i-1, S] & (roll_clear < abx.s)]
            # All column indices which roll < probability of clearance AND there is antibiotic in the previous day
            clear_idx.r <- S[abx.matrix[i-1, S] & (roll_clear < abx.r)]
            # Merge clear_idx from s abx and r abx
            clear_idx <- unique(c(clear_idx.s, clear_idx.r))
            # Clear those that pass roll and use abx to ss
            colo.matrix[i, clear_idx] <- "ss"
            # Removed those that have been cleared by abx from list of S indices
            S <- S[!(S %in% clear_idx)]
            # Roll a random number for each remaining S for chance of transmission
            roll_trans <- runif(length(S), 0, 1)
            r_idx <- S[roll_trans < prob_r]
            same_idx <- S[roll_trans >= prob_r]
            colo.matrix[i, r_idx] <- "Sr"
            colo.matrix[i, same_idx] <- "S"
        }
        
        # Update ss
        ss <- which(prev_step == "ss")
        ss <- ss[!(ss %in% already_filled)]
        if(length(ss)){
            # roll for transmission of R
            prob_r <- 1-((1-pi_r2)^r_num)
            # roll for repop of S
            prob_s <- repop.s1
            
            #print(paste(pi_r2, r_num, repop.s1, prob_r, prob_s))
            # as a Gillespie approximation probability of r and s transmission should be small enough that they do not add to 1
            if(!(prob_r+prob_s < 1)){
                stop(paste("Error stopifnot: prob_r + prob_s < 1. prob_r:", prob_r, "prob_s:", prob_s))
            }
            
            roll_r <- runif(length(ss), 0, 1)
            r_idx <- ss[roll_r < prob_r]
            roll_ss <- runif(length(ss), 0, 1)
            s_idx <- ss[(roll_ss >= prob_r) & (roll_ss < (prob_s+prob_r)) & !abx.matrix[i, ss]]
            same_idx <- ss[!(ss %in% c(r_idx, s_idx))]
            
            colo.matrix[i, s_idx] <- "S"
            colo.matrix[i, r_idx] <- "sr"
            colo.matrix[i, same_idx] <- "ss"
        }
        
        # Update Sr
        Sr <- which(prev_step == "Sr")
        Sr <- Sr[!(Sr %in% already_filled)]
        if(length(Sr)){
            
            # Roll for abx clearance
            roll_clear <- runif(length(Sr), 0, 1)
            # All column indices which roll < probability of clearance AND there is antibiotic in the previous day
            clear_idx.s <- Sr[abx.matrix[i-1, Sr] & (roll_clear < abx.s)]
            # All column indices which roll < probability of clearance AND there is antibiotic in the previous day
            clear_idx.r <- Sr[abx.matrix[i-1, Sr] & (roll_clear < abx.r)]
            # Merge clear_idx from s abx and r abx
            clear_idx <- unique(c(clear_idx.s, clear_idx.r))
            # Clear those that pass roll and use abx to ss
            colo.matrix[i, clear_idx] <- "sr"
            # Removed those that have been cleared by abx from list of S indices
            Sr <- Sr[!(Sr %in% clear_idx)]
            
            # Roll a random number for each remaining Sr for chance decolonization
            roll_decolonise <- runif(length(Sr), 0, 1)
            r_idx <- Sr[roll_decolonise < mu1]
            same_idx <- Sr[roll_decolonise >= mu1]
            colo.matrix[i, r_idx] <- "S"
            colo.matrix[i, same_idx] <- "Sr"
        }
        
        # Update sr
        sr <- which(prev_step == "sr")
        sr <- sr[!(sr %in% already_filled)]
        if(length(sr)){
            
            # Roll for repop
            roll_repop <- runif(length(sr), 0, 1)
            # If roll passes repop.r2 and has antibiotics, R grows
            sR_idx <- sr[abx.matrix[i-1, sr] & (roll_repop < repop.r2)]
            # Remove indices selected for repop.s2 event
            sr <- sr[!(sr %in% sR_idx)]
            # Instead, if roll does not pass repop.r2 event and passes repop.s2 with 0 antibiotics, S grows
            Sr_idx <- sr[!abx.matrix[i-1, sr] & (roll_repop < repop.s2)]
            # Remove indices already selected for repopulation event
            sr <- sr[!(sr %in% sR_idx)]
            
            # Roll a random number for each remaining sr for chance decolonization
            roll_decolonise <- runif(length(sr), 0, 1)
            ss_idx <- sr[roll_decolonise < mu2]
            same_idx <- sr[roll_decolonise >= mu2]
            colo.matrix[i, sR_idx] <- "sR"
            colo.matrix[i, Sr_idx] <- "Sr"
            colo.matrix[i, ss_idx] <- "ss"
            colo.matrix[i, same_idx] <- "sr"
        }
        
        # Update sR
        sR <- which(prev_step == "sR")
        sR <- sR[!(sR %in% already_filled)]
        if(length(sR)){
            # Roll for abx clearance
            roll_clear <- runif(length(sR), 0, 1)
            # if r abx used and pass roll for clearance
            clear_idx <- sR[(abx.matrix[i-1, sR] > 1)  & (roll_clear < abx.r)]
            # Remove indices selected for clearance
            sR <- sR[!(sR %in% clear_idx)]
            
            # Roll a random number for each remaining sR for chance to decolonize
            roll_decolonise <- runif(length(sR), 0, 1)
            sr_idx <- sR[roll_decolonise < depop.r]
            same_idx <- sR[roll_decolonise >= depop.r]
            colo.matrix[i, clear_idx] <- "sr"
            colo.matrix[i, sr_idx] <- "sr"
            colo.matrix[i, same_idx] <- "sR"
        }
    }
    
    return(colo.matrix)
}

diff_prevalence <- function(n.bed, mean.max.los, p.s, p.r.day1, p.r.dayafter,
                            prob_StartBact_R, prop_S_nonR, prop_Sr_inR, prop_sr_inR,
                            pi_r1, bif, mu1, mu2, abx.r, abx.s,
                            repop.r1, repop.r2, repop.s1, repop.s2, depop.r,
                            short_dur.s, long_dur.s, short_dur.r, long_dur.r, sdDur){
    
    old <- Sys.time() # get start time
    # DEBUG
    print(paste(n.bed, mean.max.los, p.s, p.r.day1, p.r.dayafter,
                prob_StartBact_R, prop_S_nonR, prop_Sr_inR, prop_sr_inR,
                pi_r1, bif, mu1, mu2, abx.r, abx.s,
                repop.r1, repop.r2, repop.s1, repop.s2, depop.r,
                short_dur.s, long_dur.s, short_dur.r, long_dur.r, sdDur))
    
    timestep <- 1
    n.day <- 500
    iterations <- 100
    iter_totalsR <- matrix(NA, nrow = n.day, ncol = iterations)

    
    for(iter in 1:iterations){
        patient.matrix <- patient.table(n.bed, n.day, mean.max.los, timestep)
        los.array <- summary.los(patient.matrix)
        abx.matrix <- abx.table(patient.matrix, los.array, p.s=p.s, p.r.day1=p.r.day1, p.r.dayafter=p.r.dayafter,
                                meanDur.s=short_dur.s, meanDur.r=short_dur.r, sdDur=sdDur, timestep=timestep)
        colo.matrix <- colo.table(patient.matrix=patient.matrix, los=los.array, 
                                  prob_StartBact_R, prop_S_nonR, prop_Sr_inR, prop_sr_inR)
        
        colo.matrix_filled_iter <- nextDay(patient.matrix, abx.matrix, colo.matrix, 
                                          pi_r1, bif, mu1, mu2, repop.r1, repop.r2, 
                                          repop.s1, repop.s2, depop.r, abx.r, abx.s, timestep)
        #Summary
        df <- data.frame(colo.matrix_filled_iter)
        iter_totalsR[, iter] <- rowSums(df == "sR")
        #print("end iteration loop")
    }
    totalsR_short <- mean(rowSums(iter_totalsR[ceiling(n.day*1/3):nrow(iter_totalsR),])/iterations/n.bed)
    
    iter_totalsR <- matrix(NA, nrow = n.day, ncol = iterations)
    for(iter in 1:iterations){
        patient.matrix <- patient.table(n.bed, n.day, mean.max.los, timestep)
        los.array <- summary.los(patient.matrix)
        abx.matrix <- abx.table(patient.matrix, los.array, p.s=p.s, p.r.day1=p.r.day1, p.r.dayafter=p.r.dayafter,
                                meanDur.s=long_dur.s, meanDur.r=long_dur.r, sdDur=sdDur, timestep=timestep)
        colo.matrix <- colo.table(patient.matrix=patient.matrix, los=los.array, 
                                  prob_StartBact_R, prop_S_nonR, prop_Sr_inR, prop_sr_inR)
        
        colo.matrix_filled_iter <- nextDay(patient.matrix, abx.matrix, colo.matrix, 
                                          pi_r1, bif, mu1, mu2, repop.r1, repop.r2, 
                                          repop.s1, repop.s2, depop.r, abx.r, abx.s, timestep)
        #Summary
        df <- data.frame(colo.matrix_filled_iter)
        iter_totalsR[,iter] <- rowSums(df == "sR")
        #print("end iteration loop")
    }
    totalsR_long <- mean(rowSums(iter_totalsR[ceiling(n.day*1/3):nrow(iter_totalsR),])/iterations/n.bed)
    
    #print(paste("totalsR_long", totalsR_long, "totalsR_short", totalsR_short))
    # print elapsed time
    new <- Sys.time() - old # calculate difference
    print(new) # print in nice format
    
    return(totalsR_long - totalsR_short)
}

parameters_binary <- c("n.bed", "mean.max.los", "p.s", "p.r.day1", "p.r.dayafter",
                      "prob_StartBact_R", "prop_S_nonR", "prop_Sr_inR", "prop_sr_inR",
                      "pi_r1", "bif", "mu1", "mu2", "abx.r", "abx.s",
                      "repop.r1", "repop.r2", "repop.s1", "repop.s2", "depop.r",
                      "short_dur.s", "long_dur.s", "short_dur.r", "long_dur.r", "sdDur")
