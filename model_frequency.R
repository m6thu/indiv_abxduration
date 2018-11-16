# 
# positive_norm_sample <- function(mean, sd){                      
#     v<-round(rnorm(1, mean=mean, sd=sd))
#     if (v < 0) {
#         v<-v*-1
#     }
#     return(v)
# }

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
    
    vector.los <- rep(0, n.days)
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

#3. Generate baseline carriage status (define function)

# Defaults from Rene's data ini_16S_log
gen_StartBact <- function(los, K, t_mean = 4.0826, t_sd = 1.1218, r_mean =1.7031, r_sd = 1.8921, n.beds, n.days){
    
    # Perhaps add check here to make sure the K is somewhat sane?
    
    number_of_patients <- dim(los)[2]
    # Unit test example: mean(t_mean), sd(t_sd)
    total_bact <- rnorm(number_of_patients, t_mean, t_sd)
    r_bact <- rnorm(number_of_patients, r_mean, r_sd)
    
    # since both are on log scale, to sum them together would need logsumexp()
    s_bact <- exp(total_bact) - exp(r_bact)
    # spin until no minus values for log... probably not best way
    while(min(s_bact) < 0){
        spin_n <- sum(s_bact < 0)
        spin_total <- rnorm(spin_n, t_mean, t_sd)
        spin_r <- rnorm(spin_n, r_mean, r_sd)
        s_bact[s_bact < 0] <- exp(spin_total) - exp(spin_r)
    }
    # convert back to log, log part of logsumexp()
    s_bact <- log(s_bact)
    
    S_Bactlevelstart <- matrix(NA, n.days, n.bed)
    R_Bactlevelstart <- matrix(NA, n.days, n.bed)
    # pad with NAs
    end_idx <- 1
    for(i in 1:number_of_patients){
        S_Bactlevelstart[end_idx:(end_idx + los[2, i] - 1)] <- c(s_bact[i], rep(NA, los[2, i]-1))
        R_Bactlevelstart[end_idx:(end_idx + los[2, i] - 1)] <- c(r_bact[i], rep(NA, los[2, i]-1))
        end_idx = end_idx + los[2, i]
    }
    
    return(list(S_Bactlevelstart, R_Bactlevelstart))
}

# 4. Update values for every day (define function)
nextDay <- function(bed_table, array_LOS, treat_table, colo_table, 
                    pi_r1, pi_r2, mu1, mu2, abx.r, abx.s,
                    repop.r1, repop.r2, repop.r3, repop.s1, repop.s2){
    
    S_table <- colo_table[[1]]
    R_table <- colo_table[[2]]
    
    # For each day (first day should be filled)
    for(i in 2:nrow(bed_table)){
        # calculate how many people has R above transmission level
        r_num <- sum(R_table[i-1,] > R_thres)
        # from number of r, calculate probability of transmission
        prob_r <- 1-((1-pi_r2)^r_num)
        # roll for transmission
        roll <- runif(1, 0, 1)
        #for each person:
        for(j in 1:ncol(bed_table)){
            if(is.na(R_table[i, j])){ # pick any; S and R should be filled in same slots
                # calculate effect of R logistic bacteria growth 
                R_grow = r_growth*R_table[i-1, j]*(1 - (R_table[i-1, j] + S_table[i-1, j])/bact_slots)
                # add effect of transmission if roll pass prob check and if previous R level is 0
                R_trans = r_trans*((roll > prob_r) & !R_table[i-1, j])
                # add effect of abx death if treat_table is r abx (== 2)
                R_abx = -(treat_table[i-1, j] > 1)*abxr_killr
                # apply effects to current table
                R_table[i, j] = R_table[i-1, j] + R_grow + R_trans + R_abx
                # trim if value goes beyond range
                if(R_table[i, j] > bact_slots){
                    R_table[i, j] = bact_slots
                }
                if(R_table[i, j] < 0){
                    R_table[i, j] = 0
                }
                
                # calculate effect of S logistic bacteria growth 
                S_grow = r_growth*S_table[i-1, j]*(1 - (R_table[i-1, j] + S_table[i-1, j])/bact_slots)
                # calculate effect of death antibiotics R and effect of death by abx S
                S_abx_s = -(treat_table[i-1, j] == 1)*abxs_kills
                S_abx_r = -(treat_table[i-1, j] > 1)*abxr_kills
                # apply effects
                S_table[i, j] = S_table[i-1, j] + R_grow + S_abx_s + S_abx_r
                # trim range
                if(S_table[i, j] > bact_slots){
                    S_table[i, j] = bact_slots
                }
                if(S_table[i, j] < 0){
                    S_table[i, j] = 0
                }
                
            }
        }
        
    }
    
    return(list(S_table, R_table))
}

diff_prevalence <- function(n.bed, n.days, mean.max.los, p.s, p.r,
                        prob_StartBact, pi_r1, pi_r2, mu1, mu2, abx.r, abx.s,
                        repop.r1, repop.r2, repop.r3, repop.s1, repop.s2,
                        iterations=10, short_dur, long_dur){
    
    iter_totalsR <- matrix(NA, nrow = n.days, ncol = iterations)
    for(iter in 1:iterations){
        
        #print(paste("iter:", iter, "y:", y_count, '-', y, "x", x_count, '-', x))
        #Generate length of stay and antibiotic duration table
        abx_iter <- abx.table(n.bed=n.bed, n.days=n.days, mean.max.los=mean.max.los, p.s=p.s, p.r=p.r, meanDur=short_dur)
        #Generate baseline carriage status
        array_LOS_iter <- array_LOS_func(los_duration=abx_iter[[2]])
        #Update values for every day
        array_StartBact_iter <- gen_StartBact(los=array_LOS_iter, prob_StartBact)
        #output
        colo_table_filled_iter <- nextDay(bed_table= abx_iter[[1]], array_LOS=array_LOS_iter, 
                                          treat_table=abx_iter[[3]], colo_table=array_StartBact_iter, 
                                          pi_r1=pi_r1, pi_r2= pi_r2, mu1=mu1, mu2=mu2, 
                                          abx.r=abx.r,abx.s=abx.s,
                                          repop.r1 = repop.r1, repop.r2 = repop.r2, repop.r3 = repop.r3, 
                                          repop.s1 = repop.s1, repop.s2 = repop.s2)
        #Summary
        df <- data.frame(colo_table_filled_iter[[2]])
        iter_totalsR[, iter] <- rowSums(df)
        #print("end iteration loop")
    }
    totalsR_short <- mean(rowSums(iter_totalsR)/iterations/n.bed)
    
    iter_totalsR <- matrix(NA, nrow = n.days, ncol = iterations)
    for(iter in 1:iterations){
        
        #print(paste("iter:", iter, "y:", y_count, '-', y, "x", x_count, '-', x))
        #Generate length of stay and antibiotic duration table
        abx_iter <- abx.table(n.bed=n.bed, n.days=n.days, mean.max.los=mean.max.los, p.s=p.s, p.r=p.r, meanDur=long_dur)
        #Generate baseline carriage status
        array_LOS_iter <- array_LOS_func(los_duration=abx_iter[[2]])
        #Update values for every day
        array_StartBact_iter <- gen_StartBact(los=array_LOS_iter, prob_StartBact)
        #output
        colo_table_filled_iter <- nextDay(bed_table= abx_iter[[1]], array_LOS=array_LOS_iter, 
                                          treat_table=abx_iter[[3]], colo_table=array_StartBact_iter, 
                                          pi_r1=pi_r1, pi_r2= pi_r2, mu1=mu1, mu2=mu2, 
                                          abx.r=abx.r,abx.s=abx.s,
                                          repop.r1 = repop.r1, repop.r2 = repop.r2, repop.r3 = repop.r3, 
                                          repop.s1 = repop.s1, repop.s2 = repop.s2)
        #Summary
        df <- data.frame(colo_table_filled_iter[[2]])
        iter_totalsR[,iter] <- rowSums(df)
        #print("end iteration loop")
    }
    totalsR_long <- mean(rowSums(iter_totalsR)/iterations/n.bed)
    
    print(paste("totalsR_long", totalsR_long, "totalsR_short", totalsR_short))
    
    return(totalsR_long - totalsR_short)
}

# whole_model(n.bed, n.days, mean.max.los, p.s, p.r,
# prob_StartBact, pi_r1, pi_r2, mu1, mu2, abx.r, abx.s,
# repop.r1, repop.r2, repop.r3, repop.s1, repop.s2,
# iterations=10, 4, 14)

######################### Heatmap Prototype (transmission vs prescription) ###########################
# 
# save_runs <- list()
# 
# ## Change parameters here: START
# y_name <- "bif"
# bif_seq <- seq(0, 500, 100)
# y_seq <- bif_seq
# 
# x_name <- "repop.s"
# repop.s_seq <- seq(0, 0.9, 0.1)
# x_seq <- repop.s_seq
# 
# iterations <- 10 # repeats per each parameter setting
# # Change parameter: END
# # Other things to consider: if any parameters are dependent add them inside loop below
# # Parameters not added here are considered fixed as defined above in Fixed Parameters
# 
# totalS <- matrix(NA, nrow = length(y_seq), ncol = length(x_seq))
# totalsr <- matrix(NA, nrow = length(y_seq), ncol = length(x_seq))
# totalSr <- matrix(NA, nrow = length(y_seq), ncol = length(x_seq))
# totalsR <- matrix(NA, nrow = length(y_seq), ncol = length(x_seq))
# 
# i <- 1 # total count
# y_count <- 1
# for(y in y_seq){
#     # calculate parameter dependencies on y paramter here
#     pi_r2 <- pi_r1 * y                  # pi_r2= probability of R transmitting to s to become sr 
#     #                                       (pi_r1 < pi_r2 if being colonised with S protects colonisation by R)
#     x_count <- 1 
#     for(x in x_seq){
#         # calculate parameter dependencies on y parameter here
#         repop.s1 <- x                        # probability of repopulation of s to become S 
#         repop.s2 <- x                          # probability of repopulation of sr to become SR 
#         repop.s3 <- x                       # probability of repopulation of sR to become sr
#         repop.r1 <- x*10                          # probability of repopulation of Sr to become sR 
#         repop.r2 <- x*10                          # probability of repopulation of sr to become sR 
#         
#         # saves for iterations
#         abx_iter <- list()
#         array_LOS_iter <- list()
#         array_StartBact_iter <- list()
#         colo_table_filled_iter <- list()
#         
#         iter_totalS <- matrix(NA, nrow = n.days, ncol = iterations)
#         iter_totalsr <- matrix(NA, nrow = n.days, ncol = iterations)
#         iter_totalSr <- matrix(NA, nrow = n.days, ncol = iterations)
#         iter_totalsR <- matrix(NA, nrow = n.days, ncol = iterations)
#         
#         for(iter in 1:iterations){
#             
#             print(paste("iter:", iter, "y:", y_count, '-', y, "x", x_count, '-', x))
#             #Generate length of stay and antibiotic duration table
#             abx_iter[[iter]] <- abx.table(n.bed, n.days, mean.max.los, p.s, p.r, meanDur=mean_dur)
#             #Generate baseline carriage status
#             array_LOS_iter[[iter]] <- array_LOS(los_duration=abx_iter[[iter]][2])
#             #Update values for every day
#             array_StartBact_iter[[iter]] <- gen_StartBact(los=array_LOS_iter[[iter]], prob_StartBact)
#             #output
#             colo_table_filled_iter[[iter]] <- nextDay(bed_table= abx_iter[[iter]][[1]], array_LOS=array_LOS_iter[[iter]], 
#                                                       treat_table=abx_iter[[iter]][[3]], colo_table=array_StartBact_iter[[iter]], 
#                                                       pi_r1=pi_r1, mu1=mu1, mu2=mu2, pi_r2=pi_r2, repop.r1 = repop.r1, 
#                                                       repop.r2 = repop.r2, repop.s1 = repop.s1, repop.s2 = repop.s2)
#             #Summary
#             df <- colo_table_filled_iter[[iter]]
#             iter_totalS[,iter] <- rowSums(df == "S")
#             iter_totalsr[,iter] <- rowSums(df == "sr")
#             iter_totalSr[,iter] <- rowSums(df == "Sr")
#             iter_totalsR[,iter] <- rowSums(df == "sR")
#             #print("end iteration loop")
#         }
#         save_runs[[i]] <- list(abx_iter, array_LOS_iter, array_StartBact_iter, colo_table_filled_iter, 
#                                iter_totalS, iter_totalsr, iter_totalSr, iter_totalsR)
#         # increment index
#         i <- i + 1
#         
#         totalS[y_count, x_count] <- mean(rowSums(iter_totalS)/iterations/n.bed)
#         totalsr[y_count, x_count] <- mean(rowSums(iter_totalsr)/iterations/n.bed)
#         totalSr[y_count, x_count] <- mean(rowSums(iter_totalSr)/iterations/n.bed)
#         totalsR[y_count, x_count] <- mean(rowSums(iter_totalsR)/iterations/n.bed)
#         
#         x_count <- x_count + 1
#         #print("end p.s loop")
#     }
#     y_count <- y_count + 1
#     #print("end pi_r1 loop")
# }
# 
# # Clean run saved to 1443_2Oct2018.RData
# image_name <- paste0(y_name, "_vs_", x_name, "_", format(Sys.time(), "%d%b%Y_%H%M%Z"))
# save.image(paste0(image_name, ".Rdata"))
# 
# require(fields)
# 
# pdf(image_name)
# par(mfrow=c(2,2))
# image.plot(t(totalS), xlab=x_name, ylab=y_name, axes=F, col=rev(terrain.colors(100)))
# title("totalS")
# axis(1, at = (1:ncol(totalS)-1)/10, labels=as.character(x_seq))
# axis(2, at = (1:nrow(totalS)-1)/10, labels=as.character(y_seq))
# image.plot(t(totalsr), xlab=x_name, ylab=y_name, axes=F, col=rev(terrain.colors(100)))
# title("totalsr")
# axis(1, at = (1:ncol(totalsr)-1)/10, labels=as.character(x_seq))
# axis(2, at = (1:nrow(totalsr)-1)/10, labels=as.character(y_seq))
# image.plot(t(totalSr), xlab=x_name, ylab=y_name, axes=F, col=rev(terrain.colors(100)))
# title("totalSr")
# axis(1, at = (1:ncol(totalSr)-1)/10, labels=as.character(x_seq))
# axis(2, at = (1:nrow(totalSr)-1)/10, labels=as.character(y_seq))
# image.plot(t(totalsR), xlab=x_name, ylab=y_name, axes=F, col=rev(terrain.colors(100)))
# title("totalsR")
# axis(1, at = (1:ncol(totalsR)-1)/10, labels=as.character(x_seq))
# axis(2, at = (1:nrow(totalsR)-1)/10, labels=as.character(y_seq))
# dev.off()
