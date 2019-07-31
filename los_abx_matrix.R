# test.abx.r.after=function(cum.r.1, all_los){
# 
#      all_admission_days= lapply(all_los, function(x) seq(3, x, by=1)) #abx r only can start after 48h
#      abx.r.after=list()
# 
#      for (i in 1:length(all_los)){
#          probs=rnorm(100) #randomly draw 100 probabilities from a normal distribution
#          probs.normalized = (probs - min(probs))/(max(probs)- min(probs)) #normalised probabilities
#          p=ecdf(probs.normalized) #cumulative distribution
#          prob.r.after=c(0,0, p((1/cum.r.1)*all_admission_days[[i]]))
#          abx.r.after.binary=rbinom(all_los[i],1,prob = prob.r.after)
# 
#          if (sum(abx.r.after.binary)>0) { #if abx.r.after is not NA i.e. abx.r started after admission
#              abx.r = which(abx.r.after.binary==1)
#              abx.r.after[[i]]= abx.r[abx.r <= all_los[i]]
#          }
#      }
#      return(abx.r.after)
#  }
#  test.abx.r.after(cum.r.1=10000, all_los=c(10,20,30,40,50,60)) 

#generate a table of number of days we want to observe (rows) -
# against number of beds in the ward (columns), filled in with patient id numbers
# accomodating for increase in length of stay with hospital acquired infections

seq_along_admissiondates=function(x){
    return(1:x)
}
seq2 <- Vectorize(seq_along_admissiondates, vectorize.args = "x")
rep2 <- Vectorize(rep.int, vectorize.args = "times")

los.abx.table <- function(n.bed, n.day, mean.max.los,
                          p.infect, p.r.day1, cum.r.1, 
                          meanDur, timestep){
    # number of slots in the patient matrix 
    # (maximum number of patients if everyone stays only for 1 day)
    n.patient.max= n.bed*n.day 
    
    #vectorise the patient id to be used for filling in the patient.matrix
    patient.id = 1:n.patient.max
    
    #length of stay for each patient if no hospital acquired infection 
    # (truncated normal distribution)
    all_los = ceiling(rexp(n.patient.max, 1/mean.max.los))
    all_los[all_los > 5*mean.max.los] = mean.max.los
    
    #decide if patients are on antibiotics 
    all_abx= rep(list(NA),n.patient.max)
    
    #day 1 of admission
    #number of days of s antibiotic is randomly drawn from a truncated normal distribution
    abx_days.s = round(rtnorm(100, mean=meanDur, lower=1))
    #number of days of r antibiotic is randomly drawn from a truncated normal distribution
    abx_days.r = round(rtnorm(100, mean=meanDur, lower=1))
    
    no.abx.r.day1=round(n.patient.max*p.infect*p.r.day1)#broad spectrum spectrum antibiotics for community acquired infection 
    no.abx.s.day1=round(n.patient.max*p.infect*(1-p.r.day1)) #narrow spectrum antibiotics for community acquired infection 
    no.abx.none.day1=n.patient.max-no.abx.r.day1-no.abx.s.day1
    id.abx=split(patient.id, sample(rep(1:3, c(no.abx.r.day1, 
                                               no.abx.s.day1,
                                               no.abx.none.day1))))
    id.abx.r.day1=id.abx$`1`
    id.abx.s.day1=id.abx$`2`
    
    all_abx[id.abx.r.day1]= rep2(2,times=sample(abx_days.r,length(id.abx.r.day1), replace = T))
    all_abx[id.abx.s.day1]= rep2(1,times=sample(abx_days.s,length(id.abx.s.day1), replace = T))
    for (i in 1:length(all_abx)){
        all_abx[[i]]=all_abx[[i]][1:all_los[i]] #extend length to los ##SLOW
        all_abx[[i]][is.na(all_abx[[i]])]=0 #those NA to be replaced by 0 
    }
    
    #risk of acquiring HAI increases with length of stay 
    # (cumulative normal distribution)
    all_admission_days = seq2(all_los) #patients' admission dates in sequence
    
    probs=rnorm(10000) #randomly draw 10000 probabilities from a normal distribution
    probs.normalized = (probs - min(probs))/(max(probs)- min(probs)) #normalised probabilities
    p=ecdf(probs.normalized) #cumulative distribution ##SLOW
    
    for (i in 1:length(all_los)){ # for every patient 
        # probs=rnorm(100) #randomly draw 100 probabilities from a normal distribution
        # probs.normalized = (probs - min(probs))/(max(probs)- min(probs)) #normalised probabilities
        # p=ecdf(probs.normalized) #cumulative distribution
        prob.r.after= p((1/cum.r.1)*all_admission_days[[i]])
        prob.r.after= c(0,0,prob.r.after [-c(1,2)])#abx r only can start after 48h
        abx.r.after.binary= rbinom(all_los[i],1, prob = prob.r.after)
        
        if (sum(abx.r.after.binary)>0) { #if abx.r.after is not NA i.e. abx.r started after admission
            abx.r.after= which(abx.r.after.binary==1) 
            abx.r.after= abx.r.after[abx.r.after <= all_los[i]]
            for (k in 1: length(abx.r.after)){
                r.dur=sample(abx_days.r,1)
                start=abx.r.after[k] #abx.r start date
                end= abx.r.after[k]+r.dur-1 #abx.r end date
                all_abx[[i]][start:end]=2
            }
        }
    }
    
    #PATIENT MATRIX
    #make up a matrix of number of days we want to observe (rows) -
    #against number of beds in the ward (columns)
    all_los= lengths(all_abx) #update length of stay incorporating abx.r
    sum_los = cumsum(all_los) #cummulative stay of patients 
    patient.matrix = matrix(NA, nrow=n.day, ncol=n.bed)
    idx = 1
    for(j in 1:n.bed){
        los_idx = suppressWarnings(max(which(sum_los < n.day))) #Suppress warning that it creates -Inf
        #creates -Inf when cummulated sum of length of stay for all patients is less than number of observational days 
        
        if(los_idx == -Inf){ # Handle case where first patient stays the whole observation duration
            #print(idx:(idx+length(los)))
            patient.matrix[, j] = rep(idx, n.day)
            #print('pat')
            #print(patient.matrix[, j])
            idx = idx+1
            all_los = all_los[-1]
        }else{
            los = all_los[1:los_idx]
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


