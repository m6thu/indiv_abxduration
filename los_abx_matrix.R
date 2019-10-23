#generate a table of number of days we want to observe (rows) -
# against number of beds in the ward (columns), filled in with patient id numbers
# accomodating for increase in length of stay with hospital acquired infections

load("norm_ecdf.Rdata")
seq_along_admissiondates <- function(x){
    return(1:x)
}
seqVec <- Vectorize(seq_along_admissiondates, vectorize.args = "x")
repVec <- Vectorize(rep.int, vectorize.args = "times")

los.abx.table <- function(n.bed, n.day, max.los,
                          p.infect, p.r.day1, cum.r.1, 
                          meanDur, timestep){
    
    # Debug
    n.bed <- 5
    n.day <- 10
    mean.max.los <- 3
    p.infect <- 0.5
    p.r.day1 <- 0.2
    cum.r.1 <- 10
    meanDur <- 1
    timestep <- 1
    
    # number of slots in the patient matrix 
    # (maximum number of patients if everyone stays only for 1 day)
    n.patient.max <- n.bed*n.day 
    
    #vectorise the patient id to be used for filling in the patient.matrix
    patient.id <- 1:n.patient.max
    
    #length of stay for each patient if no hospital acquired infection 
    # (truncated normal distribution)
    all_los = ceiling(rexp(n.patient.max, 1/max.los))
    all_los[all_los > 5*max.los] = max.los
    
    #decide if patients are on antibiotics 
    all_abx <- rep(list(NA), n.patient.max)
    
    #day 1 of admission
    #varying duration 
    #number of days of s antibiotic is randomly drawn from a truncated normal distribution
    #number of days of r antibiotic is randomly drawn from a truncated normal distribution
    # round is dependent on OS and IEEE754 means round to the even number, don't think will have much effect
    abx_days.total <- round(rtnorm(200, mean = meanDur, lower = 1))
    
    abx_days.s <- abx_days.total[1:100]
    abx_days.r <- abx_days.total[101:200]
    
    no.abx.r.day1 <- round(n.patient.max*p.infect*p.r.day1)#broad spectrum spectrum antibiotics for community acquired infection 
    no.abx.s.day1 <- round(n.patient.max*p.infect*(1 - p.r.day1)) #narrow spectrum antibiotics for community acquired infection 
    id.abx <- split(patient.id, sample(rep(1:3, c(no.abx.r.day1, 
                                               no.abx.s.day1,
                                               n.patient.max-no.abx.r.day1-no.abx.s.day1))))
    id.abx.r.day1 <- id.abx$`1`
    id.abx.s.day1 <- id.abx$`2`

    #uniform durations 
    all_abx[id.abx.r.day1]= rep(list(rep(2, meanDur)), length(id.abx.r.day1))
    all_abx[id.abx.s.day1]= rep(list(rep(1, meanDur)), length(id.abx.s.day1))
    
    for (i in 1:length(all_abx)){
        all_abx[[i]] <- all_abx[[i]][1:all_los[i]] #extend length to los
        all_abx[[i]][is.na(all_abx[[i]])] <- 0 #those NA to be replaced by 0 
    }
    
    #risk of acquiring HAI increases with length of stay 
    # (cumulative normal distribution)
    all_admission_days <- seqVec(all_los) #1:los vectorized, patients' admission dates in sequence
    
    #probs <- rnorm(10) #randomly draw 100 probabilities from a normal distribution
    #probs.normalized <- (probs - min(probs))/(max(probs) - min(probs)) #shift distribution to start at 0 and have 
    #p <- ecdf(probs.normalized) #cumulative distribution
    # The above method consumes time by having to be run every simulation
    # It is also sensitive to the number of rnorm drawn therefore
    # I think use the below function instead, see more notes in norm_ecdf.R -Fai
    p <- norm_ecdf # imported normalized function
    
    for (i in 1:length(all_los)){ # for every patient 
        # probs=rnorm(100) #randomly draw 100 probabilities from a normal distribution
        # probs.normalized = (probs - min(probs))/(max(probs)- min(probs)) #normalised probabilities
        # p=ecdf(probs.normalized) #cumulative distribution
        prob.r.after <- p((1/cum.r.1)*all_admission_days[[i]]) # get the cummulative probability for each day
        prob.r.after <- c(0, 0, prob.r.after[-c(1,2)]) #abx r only can start after 48h, delete first 2 cases
        abx.r.after.binary <- rbinom(all_los[i], 1, prob = prob.r.after)
        
        if (sum(abx.r.after.binary) > 0) { #if abx.r.after is not NA i.e. abx.r started after admission
            abx.r.after <- which(abx.r.after.binary==1) 
            abx.r.after <- abx.r.after[abx.r.after <= all_los[i]]
            for (k in 1: length(abx.r.after)){
                r.dur <- sample(abx_days.r,1)
                start <- abx.r.after[k] #abx.r start date
                end <- abx.r.after[k]+r.dur-1 #abx.r end date
                all_abx[[i]][start:end] <- 2
            }
        }
    }
    
    #PATIENT MATRIX
    #make up a matrix of number of days we want to observe (rows) -
    #against number of beds in the ward (columns)
    all_los <- lengths(all_abx) #update length of stay incorporating abx.r
    sum_los <- cumsum(all_los) #cummulative stay of patients 
    patient.matrix <- matrix(NA, nrow=n.day, ncol=n.bed)
    idx <- 1
    for(j in 1:n.bed){
        los_idx <- suppressWarnings(max(which(sum_los < n.day))) #Suppress warning that it creates -Inf
        #creates -Inf when cummulated sum of length of stay for all patients is less than number of observational days 
        
        if(los_idx == -Inf){ # Handle case where first patient stays the whole observation duration
            #print(idx:(idx+length(los)))
            patient.matrix[, j] <- rep(idx, n.day)
            #print('pat')
            #print(patient.matrix[, j])
            idx <- idx+1
            all_los <- all_los[-1]
        }else{
            los <- all_los[1:los_idx]
            #print(idx:(idx+length(los)))
            patient.matrix[, j] = rep(idx:(idx+length(los)), c(los, n.day-sum(los)))
            #print('pat')
            #print(patient.matrix[, j])
            idx = idx+los_idx+1
            all_los = all_los[-(1:(los_idx+1))]
        }
        sum_los = cumsum(all_los)
    }
    patient.matrix=patient.matrix[rep(1:nrow(patient.matrix), each= timestep), ]
    
    #ANTIBIOTIC MATRIX 
    #make up a matrix of number of days we want to observe (rows) -
    #against number of beds in the ward (columns)
    abx.matrix = matrix(NA, nrow=n.day, ncol=n.bed)
    sum_abx= cumsum(lengths(all_abx))
    idx = 1
    for(j in 1:n.bed){
        los_idx = suppressWarnings(max(which(sum_abx < n.day))) #the patient id that fill each bed (column) 
        #Suppress warning that it creates -Inf
        #creates -Inf when the first patient is all_los > number of observational days 
        
        if(los_idx == -Inf){ # case where first patient stays the whole observation duration
            abx.matrix[, j] = all_abx[[1]][1:n.day]
            idx = idx+1
            all_abx = all_abx[-1]
        }else{
            abx = unlist(all_abx[1:(los_idx+1)])
            abx.matrix[, j] = abx[1:n.day]
            idx = idx+los_idx+1
            all_abx = all_abx[-(1:(los_idx+1))]
        }
        sum_abx = cumsum(lengths(all_abx))
    }
    abx.matrix=abx.matrix[rep(1:nrow(abx.matrix), each = timestep), ]
    
    return(list(patient.matrix, abx.matrix))
}

#frequency summary of patient.matrix - patient id against number of days of stay for each patient
summary.los <- function(patient.matrix){    
    
    # Summarize how often each patient.id is encountered to get days for each id
    los.dur = table(patient.matrix) ##SLOW
    los_duration = array(dim = c(2, max(patient.matrix)))
    # Attach patient ID on 1st row
    los_duration[1,] = 1:max(patient.matrix)
    # Put summary of days on 2nd row
    los_duration[2,] = los.dur
    
    return(los_duration)
}


