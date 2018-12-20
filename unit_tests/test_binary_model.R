# Regression testing for binary model
rm(list=ls()) # Clean working environment
source("model_binary.R") # Load model for testing

# Please add as appropriate, be as pendantic as much as possible
# Cleanest way to run. Source above, select only block of test case, try running

# Generates an artificial patient table where each patient stays an equal number of days "each"
# matrix(rep(rep(1:100, rep(3, 100)), 20) + rep((0:19)*100, rep(300, 20)), ncol=20, nrow=300)
patient.table.equal <- function(n.bed, n.day, each){
    # cases where n.day is not divisible by each is not tested
    stopifnot(n.day %% each == 0)
    vec <- rep(rep(1:(n.day/each), rep(each, n.day/each)), n.bed) + rep((0:(n.bed-1))*(n.day/each), rep(n.day, n.bed))
    mat <- matrix(vec, ncol=n.bed, nrow=n.day)
    return(mat)
}

# Generates and artificial colonization starting table that uses probability as fixed proportion instead
# meant to be used with patient.table.equal(... , each=2)
colo.table.prop <- function(patient.matrix, los.array, prob_StartBact_R, prop_S_nonR){
    
    each <- ncol(patient.matrix)
    r_num <- round(each*prob_StartBact_R)
    s_num <- round(each*prop_S_nonR*(1-prob_StartBact_R))
    ss_num <- each - r_num - s_num
    vec <- c(rep("R", r_num), rep("S", s_num), rep("ss", ss_num))
    colo.matrix <- t(matrix(rep(c(vec, rep(NA, each)), nrow(patient.matrix)/2), 
                            ncol=nrow(patient.matrix), nrow=ncol(patient.matrix)))
    return(colo.matrix)
}

# Template
# ---------------------------------------------------------------------------------------
# Test name: 
# Test summary: 
# Function tested: 
# Pre-conditions: 
# Expected result: 

############################################# Function tests ##################################################
#################################### patient.table ####################################
# ---------------------------------------------------------------------------------------
# Test name: simple run
# Test summary: Tests patient.table functions exits with no errors, gives correct mean, and gives correct shape
# Function tested: patient.table
# Pre-conditions:
rm(list=ls()) # Clean working environment
source("model_binary.R") # Load model for testing
test_mean <- 5
tolerance <- 1
patient_mat <- patient.table(n.bed = 20, n.day = 300, mean.max.los = test_mean, timestep=1)
# Expected result: at high points, patient table should give exponential distribution with mean given within tolerance
hist(table(patient_mat)) # eyeball that this looks like an exponential distribution and that maximum value at tail makes sense
# Expected result: mean from table is the same mean given within tolerance
stopifnot(abs(mean(table(patient_mat)) - test_mean) < tolerance)
# Expected result: dimensions are correct dimensions
stopifnot(dim(patient_mat) == c(300, 20))

# ---------------------------------------------------------------------------------------
# Test name: simple run
# Test summary: For >1 timesteps, tests patient.table functions exits with no errors, gives correct mean, and gives correct shape
# Function tested: patient.table
# Pre-conditions:
rm(list=ls()) # Clean working environment
source("model_binary.R") # Load model for testing
test_mean <- 7
tolerance <- 1*3
patient_mat.m <- patient.table(n.bed = 100, n.day = 90, mean.max.los = test_mean, timestep=3)
# Expected result: looks like an exponential distribution and that maximum value at tail makes sense
hist(table(patient_mat.m))
# Expected result: mean from table is the same mean given within tolerance
stopifnot(abs(mean(table(patient_mat.m)) - test_mean*3) < tolerance)
# Expected result: dimensions are correct dimensions
stopifnot(dim(patient_mat.m) == c(270, 100))

#################################### summary.los ####################################
# ---------------------------------------------------------------------------------------
# Test name: 
# Test summary: 
# Function tested: summary.los
# Pre-conditions: 
rm(list=ls()) # Clean working environment
source("model_binary.R") # Load model for testing
test_mean <- 3
tolerance <- 1
patient_mat <- patient.table(n.bed = 20, n.day = 300, mean.max.los = test_mean, timestep=1)
los_duration <- summary.los(patient_mat)
# Expected result: 
hist(los_duration[2, ]) # check distribution of length of stay, should be exponential
# Expected result: 
hist(table(patient_mat)) # chould be the same as input, los.array should = table(patient.table.output)
# Expected result: 
stopifnot(max(patient_mat.s) == max(los_duration)) # highest number from array should be exactly the same as highest id in patient.matrix

# ---------------------------------------------------------------------------------------
# Test name: 
# Test summary: 
# Function tested: summary.los
# Pre-conditions: 
rm(list=ls()) # Clean working environment
source("model_binary.R") # Load model for testing
test_mean <- 3
tolerance <- 1*3
patient_mat.m <- patient.table(n.bed = 100, n.day = 90, mean.max.los = test_mean, timestep=3)
los_duration.m <- summary.los(patient_mat.m)
# Expected result: 
hist(los_duration.m[2, ]) # check distribution of length of stay, should be exponential
# Expected result: 
hist(table(patient_mat.m)) # chould be the same as input, los.array should = table(patient.table.output)
# Expected result: 
stopifnot(max(patient_mat.m) == max(los_duration.m)) # highest number from array should be exactly the same as highest id in patient.matrix

#################################### Test abx matrix generation ####################################
# ---------------------------------------------------------------------------------------
# Test name: simple run
# Test summary: equivalent to simple_model test
# Function tested: abx.table
# Pre-conditions:
rm(list=ls()) # Clean working environment
source("model_binary.R") # Load model for testing
test_p <- 0.3
test_mean <- 10
tolerance <- 1
patient_mat.s <- patient.table(n.bed = 100, n.day = 90, mean.max.los = test_mean, timestep=1)
los_duration.s <- summary.los(patient_mat.s)
# Expected result: Should be equivalent to simple_model test
abx.matrix.s <- abx.table(patient.matrix=patient_mat.s, los.array=los_duration.s, p.s=test_p, p.r.day1=0, p.r.dayafter=0.1,
                          meanDur.s=test_mean, meanDur.r=test_mean, sdDur=1, timestep=1)
# To use commented test below, need to edit code: c(rep(1, abx_person), rep(0, max_days-abx_person)) to c(rep(1, abx_person), rep(1, max_days-abx_person))
# stopifnot(abs(sum(abx.matrix.s > 0)/length(abx.matrix.s) - test_p) < tolerance) # overall number of 1s approx. probability

# ---------------------------------------------------------------------------------------
# Test name: meanDur.s, p.s, matrix shape
# Test summary: 100% effectiveness, no clipping of antibiotics
# Function tested: abx.table
# Pre-conditions: 
rm(list=ls()) # Clean working environment
source("model_binary.R") # Load model for testing
test_p <- 1
test_mean <- 10
tolerance <- 1
patient_mat.s <- patient.table(n.bed = 10, n.day = 3000, mean.max.los = 100, timestep=1) # have days long enough to not create clipping
los_duration.s <- summary.los(patient_mat.s)
abx.matrix.s <- abx.table(patient_mat.s, los_duration.s, p.s=test_p, p.r.day1=0, p.r.dayafter=0,
                          meanDur.s=test_mean, meanDur.r=0, sdDur=1, timestep=1)
abx_summary <- table(abx.matrix.s*patient_mat.s)
# Expected result: 
hist(abx_summary[2:length(abx_summary)], breaks=20) # distribution of abx duration shape should be truncated norm as used
stopifnot(abs(mean(abx_summary[2:length(abx_summary)]) - test_mean) < tolerance) # mean of abx duration is within tolerance
stopifnot(dim(abx.matrix.s) == dim(patient_mat.s)) # dimensions must equal patient.matrix

# ---------------------------------------------------------------------------------------
# Test name: p.s, test_los_mean
# Test summary: # 100% effectiveness, all clipping of antibiotics
# Function tested: abx.table
# Pre-conditions: 
rm(list=ls()) # Clean working environment
source("model_binary.R") # Load model for testing
# should track exponential function of length of stay
test_p <- 1
test_mean <- 10
test_los_mean <- 5 
tolerance <- 1
patient_mat.s <- patient.table(n.bed = 20, n.day = 3000, mean.max.los = test_los_mean, timestep=1) # have days short to create clipping
los_duration.s <- summary.los(patient_mat.s)
abx.matrix.s <- abx.table(patient_mat.s, los_duration.s, p.s=test_p, p.r.day1=0, p.r.dayafter=0,
                          meanDur.s=test_mean, meanDur.r=0, sdDur=1, timestep=1)
abx_summary <- table(abx.matrix.s*patient_mat.s)
# Expected output: 
hist(abx_summary[2:length(abx_summary)], breaks=20) # distribution of abx duration shape should be exponential function, small bump at test_mean for "lucky hits"
stopifnot(abs(mean(abx_summary[2:length(abx_summary)]) - test_los_mean) < tolerance) # mean of abx duration is mean of los instead
stopifnot(dim(abx.matrix.s) == dim(patient_mat.s)) # dimensions must equal patient.matrix

# ---------------------------------------------------------------------------------------
# Test name: p.r.day1
# Test summary: test starting probability of R antibiotics at 100% effectiveness, all clipping of antibiotics, 
# Function tested: abx.table
# Pre-conditions: 
rm(list=ls()) # Clean working environment
source("model_binary.R") # Load model for testing
test_pr <- 0.6
test_mean <- 10
tolerance <- 1
patient_mat.s <- patient.table(n.bed = 10, n.day = 3000, mean.max.los = 100, timestep=1) # have days long enough to not create clipping
los_duration.s <- summary.los(patient_mat.s)
abx.matrix.s <- abx.table(patient_mat.s, los_duration.s, p.s=0, p.r.day1=test_pr, p.r.dayafter=0,
                          meanDur.s=0, meanDur.r=test_mean, sdDur=1, timestep=1)
abx_summary <- table(abx.matrix.s*patient_mat.s)
# Expected output: hould track exponential function of length of stay
hist(abx_summary[2:length(abx_summary)], breaks=20) # distribution of abx duration shape should be truncated norm as used
# Expected output: 
stopifnot(abs(mean(abx_summary[2:length(abx_summary)]) - test_mean) < tolerance) # mean of abx duration is within tolerance
# Expected output: 
stopifnot(dim(abx.matrix.s) == dim(patient_mat.s)) # dimensions must equal patient.matrix



# ---------------------------------------------------------------------------------------
# Test name: multi timestep, simple run
# Test summary: function exits cleanly with no errors
# Function tested: abx.table
# Pre-conditions: 
rm(list=ls()) # Clean working environment
source("model_binary.R") # Load model for testing
test_p <- 0.3
test_mean <- 10
tolerance <- 1
patient_mat.s <- patient.table(n.bed = 20, n.day = 3000, mean.max.los = 100, timestep=2) # have days long enough to not create clipping
los_duration.s <- summary.los(patient_mat.s)
# Expected output:
abx.matrix.s <- abx.table(patient_mat.s, los_duration.s, p.s=test_p, p.r.day1=0, p.r.dayafter=0,
                          meanDur.s=test_mean, meanDur.r=0, sdDur=1, timestep=2)
# To use commented test below, need to edit code: c(rep(1, abx_person), rep(0, max_days-abx_person)) to c(rep(1, abx_person), rep(1, max_days-abx_person))
# stopifnot(abs(sum(abx.matrix.s > 0)/length(abx.matrix.s) - test_p) < tolerance) # overall number of 1s approx. probability

# ---------------------------------------------------------------------------------------
# Test name: multi timestep p.s
# Test summary: # 100% effectiveness, all clipping of antibiotics
# Function tested: abx.table
# Pre-conditions: 
rm(list=ls()) # Clean working environment
source("model_binary.R") # Load model for testing
test_p <- 1
test_mean <- 10
tolerance <- 1*2
patient_mat.s <- patient.table(n.bed = 20, n.day = 3000, mean.max.los = 100, timestep=2) # have days long enough to not create clipping
los_duration.s <- summary.los(patient_mat.s)
abx.matrix.s <- abx.table(patient_mat.s, los_duration.s, p.s=test_p, p.r.day1=0, p.r.dayafter=0,
                          meanDur.s=test_mean, meanDur.r=0, sdDur=1, timestep=2)
abx_summary <- table(abx.matrix.s*patient_mat.s)
# Expected output: should track exponential function of length of stay
hist(abx_summary[2:length(abx_summary)], breaks=20) # distribution of abx duration shape should be truncated norm as used
stopifnot(abs(mean(abx_summary[2:length(abx_summary)]) - test_mean*2) < tolerance) # mean of abx duration is within tolerance
stopifnot(dim(abx.matrix.s) == dim(patient_mat.s)) # dimensions must equal patient.matrix

# ---------------------------------------------------------------------------------------
# Test name: 
# Test summary: # 100% effectiveness, all clipping of antibiotics
# Function tested: 
# Pre-conditions: 
# should track exponential function of length of stay
# 100% effectiveness, all clipping of antibiotics
# should track exponential function of length of stay
rm(list=ls()) # Clean working environment
source("model_binary.R") # Load model for testing
test_p <- 1
test_mean <- 10
test_los_mean <- 5 
tolerance <- 1*2
patient_mat.s <- patient.table(n.bed = 20, n.day = 3000, mean.max.los = test_los_mean, timestep=2) # have days short to create clipping
los_duration.s <- summary.los(patient_mat.s)
abx.matrix.s <- abx.table(patient_mat.s, los_duration.s, p.s=test_p, p.r.day1=0, p.r.dayafter=0,
                          meanDur.s=test_mean, meanDur.r=0, sdDur=1, timestep=2)
abx_summary <- table(abx.matrix.s*patient_mat.s)
# Expected output: 
hist(abx_summary[2:length(abx_summary)], breaks=20) # distribution of abx duration shape should be exponential function, small bump at test_mean for "lucky hits"
stopifnot(abs(mean(abx_summary[2:length(abx_summary)]) - test_los_mean*2) < tolerance) # mean of abx duration is mean of los instead
stopifnot(dim(abx.matrix.s) == dim(patient_mat.s)) # dimensions must equal patient.matrix

# ---------------------------------------------------------------------------------------
# Test name: 
# Test summary: # 100% effectiveness, all clipping of antibiotics
# Function tested: 
# Pre-conditions: 
# should track exponential function of length of stay
# Case: Multiple timestep, multiple abx prob

# ---------------------------------------------------------------------------------------
# Test name: 
# Test summary: # 100% effectiveness, all clipping of antibiotics
# Function tested: 
# Pre-conditions: 
# should track exponential function of length of stay
# Case: Single timestep, check prob of starting abx per day

# ---------------------------------------------------------------------------------------
# Test name: 
# Test summary: # 100% effectiveness, all clipping of antibiotics
# Function tested: 
# Pre-conditions: 
# should track exponential function of length of stay
# Case: Multiple timestep, check prob of starting abx per day

#################################### Test starting bacteria generation ####################################

# ---------------------------------------------------------------------------------------
# Test name: Single time step
# Test summary: # 100% effectiveness, all clipping of antibiotics
# Function tested: 
# Pre-conditions: 
# should track exponential function of length of stay
rm(list=ls()) # Clean working environment
source("model_binary.R") # Load model for testing
#prob_StartBact_R, prop_S_nonR, prop_Sr_inR, prop_sr_inR
prob_StartBact_R <- 0.2
prop_S_nonR <- 0.3
prop_Sr_inR <- 0.5
prop_sr_inR <- 0.5
tolerance <- 0.1
prob_start_S <- prop_S_nonR*(1-prob_StartBact_R)
prob_start_Sr <- prop_Sr_inR*prob_StartBact_R
prob_start_sr <- prop_sr_inR*prob_StartBact_R
patient.matrix <- patient.table(n.bed = 20, n.day = 300, mean.max.los = 5, timestep=1)
los.array <- summary.los(patient.matrix)
colo.matrix <- colo.table(patient.matrix, los.array, prob_StartBact_R, prop_S_nonR, prop_Sr_inR, prop_sr_inR)
# Expected output: 
stopifnot(abs(sum(colo.matrix == "S", na.rm=TRUE)/sum(!is.na(colo.matrix)) - prop_S_nonR*(1-prob_StartBact_R)) < tolerance) # generated R probs is correct
stopifnot(abs(sum(colo.matrix == "ss", na.rm=TRUE)/sum(!is.na(colo.matrix)) - (1-prob_start_S-prob_StartBact_R)) < tolerance) # generated R probs is correct
stopifnot(abs(sum(colo.matrix == "Sr", na.rm=TRUE)/sum(!is.na(colo.matrix)) - (prop_Sr_inR*prob_StartBact_R)) < tolerance) # generated R probs is correct
stopifnot(abs(sum(colo.matrix == "sr", na.rm=TRUE)/sum(!is.na(colo.matrix)) - (prop_sr_inR*prob_StartBact_R)) < tolerance) # generated R probs is correct
stopifnot(abs(sum(colo.matrix == "sR", na.rm=TRUE)/sum(!is.na(colo.matrix)) - (prob_StartBact_R-prob_start_Sr-prob_start_sr)) < tolerance) # generated R probs is correct
stopifnot(sum(is.na(colo.matrix[1, ])) == 0) # Top row is fully initialized
stopifnot(dim(colo.matrix) == dim(patient.matrix)) # dimension must equal patient.matrix
patient_idx <- cumsum(c(1, los.array[2,]))
stopifnot(sum(!(which(!is.na(colo.matrix)) == patient_idx[-length(patient_idx)])) == 0) # starts on same position as patient id

# ---------------------------------------------------------------------------------------
# Test name: Multiple time step
# Test summary: # 100% effectiveness, all clipping of antibiotics
# Function tested: 
# Pre-conditions: 
# should track exponential function of length of stay
rm(list=ls()) # Clean working environment
source("model_binary.R") # Load model for testing
prob_StartBact_R <- 0.2
prop_S_nonR <- 0.3
prop_Sr_inR <- 0.5
prop_sr_inR <- 0.5
tolerance <- 0.1
prob_start_S <- prop_S_nonR*(1-prob_StartBact_R)
prob_start_Sr <- prop_Sr_inR*prob_StartBact_R
prob_start_sr <- prop_sr_inR*prob_StartBact_R
patient.matrix <- patient.table(n.bed = 20, n.day = 300, mean.max.los = 5, timestep=3)
los.array <- summary.los(patient.matrix)
colo.matrix <- colo.table(patient.matrix, los.array, prob_StartBact_R, prop_S_nonR, prop_Sr_inR, prop_sr_inR)
# Expected result: 
stopifnot(abs(sum(colo.matrix == "S", na.rm=TRUE)/sum(!is.na(colo.matrix)) - prop_S_nonR*(1-prob_StartBact_R)) < tolerance) # generated R probs is correct
stopifnot(abs(sum(colo.matrix == "ss", na.rm=TRUE)/sum(!is.na(colo.matrix)) - (1-prob_start_S-prob_StartBact_R)) < tolerance) # generated R probs is correct
stopifnot(abs(sum(colo.matrix == "Sr", na.rm=TRUE)/sum(!is.na(colo.matrix)) - (prop_Sr_inR*prob_StartBact_R)) < tolerance) # generated R probs is correct
stopifnot(abs(sum(colo.matrix == "sr", na.rm=TRUE)/sum(!is.na(colo.matrix)) - (prop_sr_inR*prob_StartBact_R)) < tolerance) # generated R probs is correct
stopifnot(abs(sum(colo.matrix == "sR", na.rm=TRUE)/sum(!is.na(colo.matrix)) - (prob_StartBact_R-prob_start_Sr-prob_start_sr)) < tolerance) # generated R probs is correct
stopifnot(sum(is.na(colo.matrix[1, ])) == 0) # Top row is fully initialized
stopifnot(dim(colo.matrix) == dim(patient.matrix)) # dimension must equal patient.matrix
patient_idx <- cumsum(c(1, los.array[2,]))
stopifnot(sum(!(which(!is.na(colo.matrix)) == patient_idx[-length(patient_idx)])) == 0) # starts on same position as patient id

##################################### Test daily updates ####################################
# ---------------------------------------------------------------------------------------
# Test name: 
# Test summary: Feed only certain states, check probability of conversion to another state without abx and with abx 
# Pre-conditions: 
# Expected result: 
# ---------------------------------------------------------------------------------------
# Test name: 
# Test summary: 
# Function tested: 
# Pre-conditions: 
# Expected result: 

# nextDay <- function(patient.matrix, abx.matrix, colo.matrix, 
#                     pi_r1, bif, mu1, mu2, repop.r1, repop.r2, 
#                     repop.s1, repop.s2,repop.s3, abx.r, abx.s){
#     
#     pi_r2 <- pi_r1 * bif                 # pi_r2= probability of R transmitting to s to become sr 
#     
#     # For each day (first day should be filled)
#     for(i in 2:nrow(patient.matrix)){
#         # For each bed
#         for(j in 1:ncol(patient.matrix)){
#             #case S
#             #print(paste("i:", i, "j:", j))
#             #print(colo.matrix[i-1, j])
#             if(is.na(colo.matrix[i, j])){
#                 r_num <- 5 #sum(colo.matrix[i-1,] == "sR") #only R can be transmitted 
#                 if(colo.matrix[i-1, j] == "S"){
#                     #print("----case S")
#                     # check antibiotic
#                     roll_clear <- runif(1, 0, 1)
#                     roll_transmit <- runif(1, 0, 1)
#                     prob_r <- 1-((1-pi_r1)^r_num)
#                     if (abx.matrix[i-1, j] == 1 & roll_clear < abx.s){
#                         colo.matrix[i, j] <- "ss"
#                     } else if (abx.matrix[i-1, j] > 1 & roll_clear < abx.r){
#                         colo.matrix[i, j] <- "ss"
#                     } else if (roll_transmit < prob_r){ 
#                         colo.matrix[i, j] <- "Sr"
#                     }else {
#                         colo.matrix[i, j] <- "S"
#                     }
#                     
#                     # case ss
#                 }else if(colo.matrix[i-1, j] == "ss"){
#                     #print("----case s")
#                     # roll for transmission of r
#                     roll_r <- runif(1, 0, 1)
#                     prob_r <- 1-((1-pi_r2)^r_num)
#                     # roll for repopulation of s to become S
#                     roll_ss <- runif(1, 0, 1)
#                     if (roll_r < prob_r) { 
#                         colo.matrix[i,j]<-"sr" 
#                     } else if ( roll_ss < repop.s1) {
#                         colo.matrix[i,j]<-"S" 
#                     } else{
#                         colo.matrix[i, j] <- "ss"
#                     }
#                     
#                     # case Sr
#                 }else if(colo.matrix[i-1, j] == "Sr"){
#                     #print("----case Sr")
#                     # check antibiotics 
#                     roll_clear <- runif(1, 0, 1)
#                     roll_decolonise <- runif(1, 0, 1)
#                     if(abx.matrix[i-1, j] ==1 & roll_clear < abx.s){
#                         colo.matrix[i, j] <- "sr"
#                     } else if (abx.matrix[i-1, j] >1 & roll_clear < abx.r ) {
#                         colo.matrix[i, j] <- "sr"
#                     } else if(roll_decolonise < mu1){ 
#                         colo.matrix[i, j] <- "S"
#                     }else {
#                         colo.matrix[i, j] <- "Sr"
#                     }
#                     
#                     # case sr
#                 }else if(colo.matrix[i-1, j] == "sr"){
#                     #print("----case sr")
#                     roll_repop <- runif(1, 0, 1)
#                     roll_decolonise <- runif(1, 0, 1)
#                     # check antibiotics
#                     if(abx.matrix[i-1, j] == 1 & roll_repop < repop.r2){
#                         colo.matrix[i, j] <- "sR"
#                     }else if(abx.matrix[i-1, j] == 0 & roll_repop < repop.s2){
#                         colo.matrix[i, j] <- "Sr"
#                         # }else if(abx.matrix[i-1, j] == 0 & roll_repop < repop.r3){
#                         #     colo.matrix[i, j] <- "sR"
#                     }else if (roll_decolonise < mu2){
#                         colo.matrix[i, j] <- "ss"
#                     }else {
#                         colo.matrix[i, j] <- "sr"
#                     }
#                     
#                     # case sR
#                 }else if(colo.matrix[i-1, j] == "sR"){
#                     #print("----case sR")
#                     roll_clear <- runif(1, 0, 1)
#                     if(abx.matrix[i-1, j] > 1 & roll_clear < abx.r){
#                         colo.matrix[i, j] <- "sr"
#                     }else if (abx.matrix[i-1, j] == 0 & roll_clear < repop.s3){
#                         colo.matrix[i, j] <- "sr"
#                     }else {
#                         colo.matrix[i, j] <- "sR"
#                     }
#                 }else{
#                     print("error")
#                     colo.matrix[i, j] <- "E"
#                 }
#             } # if 0
#         }  # for j
#     } # for i
#     
#     return(colo.matrix)
# }
# test.abx.table<-abx.table(patient.matrix=test.patient.table, los.array=test.los.array, p=0.5, meanDur=4, sdDur=2, timestep=1)
# test.table<-nextDay(patient.matrix=test.patient.table, abx.matrix=test.abx.table, colo.matrix=test.colo.table, 
#                     pi_r1=0.5, bif=0, mu1=0, mu2=0, repop.r1=0, repop.r2=0, 
#                     repop.s1=0, repop.s2=0,repop.s3=0, abx.r=0, abx.s=0)
# #check pi_r1
# check_pi_r1<- 0
# for (i in 2:nrow(test.table)) {
#     for (j in 1:ncol(test.table)) {
#         if (test.table[i,j]=="Sr" & test.table[i-1,j]=="S" & is.na(test.colo.table[i,j])) {
#             check_pi_r1<- check_pi_r1+1
#         }
#     }
# }
# check_abx<-0
# for (i in 2:nrow(test.table)) {
#     for (j in 1:ncol(test.table)) {
#         if (test.table[i,j]=="ss" & test.table[i-1,j]=="S" & is.na(test.colo.table[i,j])) {
#             check_abx<- check_abx+1
#         }
#     }
# }
# check_S<-0
# for (i in 2:nrow(test.table)) {
#     for (j in 1:ncol(test.table)) {
#         if (test.table[i,j]=="S" & test.table[i-1,j]=="S" & is.na(test.colo.table[i,j])) {
#             check_S<- check_S+1
#         }
#     }
# }
# test.pi_r1<-check_pi_r1/sum(check_abx,check_S,check_pi_r1) #get pi_r1
# 1-((1-test.pi_r1)^1/5)

############################################ Integration tests ##################################################
# Test diff_prevalence
diff_prevalence(n.bed=20, mean.max.los=4, p.s=0.1, p.r.day1=0.2, p.r.dayafter=0.01,
                prob_StartBact_R=0.3, prop_S_nonR=0.1, prop_Sr_inR=0.1, prop_sr_inR=0.1,
                pi_r1=0.1, bif=2, mu1=0.10, mu2=0.10, abx.r=0.10, abx.s=0.10,
                repop.r1=0.10, repop.r2=0.10, repop.s1=0.10, repop.s2=0.10, depop.r=0.10,
                short_dur.s=2, long_dur.s=10, short_dur.r=2, long_dur.r=10, sdDur=5)


# Case: At no transmission and no abx
# Expected output: base level of starting Rs generated

# Case: No transmission with abx
# Expected output:

# Case: Transmission with no abx
# Expected output: 


# Test scenario

# Case: Abx on everyone
# Expected output: No S at all

# Case: High transmission
# Expected output: ?

############################################# Debug rigs ################################################
# Uncomment when needed (so they don't interfere with regression tests)
# ---------------------------------------------------------------------------------------
# Test name: Rig for closing issue #24 on git
# Test summary: Bug fixed
# Function tested: abx.table
# Pre-conditions:
# rm(list=ls()) # Clean working environment
# source("model_binary.R") # Load model for testing
# diff_prevalence(6.6875, 8.1625, 0.016125, 0.469375, 0.431251375,
#                 0.14475, 0.04812875, 0.012375, 0.019875,
#                 0.332625, 0.332625, 0.555125, 8.75005625, 0.022125, 0.004875,
#                 0.020625, 0.6625, 0.83075, 0.12025, 0.46325,
#                 4.15, 5.95, 10.9625, 20.8625, 1.15)
# 
# # n.bed, mean.max.los, p.s, p.r.day1, p.r.dayafter,
# # prob_StartBact_R, prop_S_nonR, prop_Sr_inR, prop_sr_inR,
# # pi_r1, bif, mu1, mu2, abx.r, abx.s,
# # repop.r1, repop.r2, repop.s1, repop.s2, depop.r,
# # short_dur.s, long_dur.s, short_dur.r, long_dur.r, sdDur
# n.bed=6.6875; mean.max.los=8.1625; p.s=0.016125; p.r.day1=0.469375; p.r.dayafter=0.431251375;
# prob_StartBact_R=0.14475; prop_S_nonR=0.04812875; prop_Sr_inR=0.012375; prop_sr_inR=0.019875; 
# pi_r1=0.332625; bif=0.332625; mu1=0.555125; mu2=8.75005625; abx.r=0.022125; abx.s=0.004875;
# repop.r1=0.020625; repop.r2=0.6625; repop.s1=0.83075; repop.s2=0.12025; depop.r=0.46325;
# short_dur.s=4.15; long_dur.s=5.95; short_dur.r=10.9625; long_dur.r=20.8625; sdDur=1.15
# timestep <- 1
# n.day <- 500
# iterations <- 100
# patient.matrix <- patient.table(n.bed=n.bed, n.day, mean.max.los=mean.max.los, timestep)
# los.array <- summary.los(patient.matrix)
# abx.matrix <- abx.table(patient.matrix, los.array, p.s=p.s, p.r.day1=p.r.day1, p.r.dayafter=p.r.dayafter,
#                         meanDur.s=short_dur.s, meanDur.r=short_dur.r, sdDur=sdDur, timestep=timestep)
# colo.matrix <- colo.table(patient.matrix=patient.matrix, los=los.array, 
#                           prob_StartBact_R, prop_S_nonR, prop_Sr_inR, prop_sr_inR)
# 
# colo.matrix_filled_iter <- nextDay(patient.matrix, abx.matrix, colo.matrix, 
#                                    pi_r1, bif, mu1, mu2, repop.r1, repop.r2, 
#                                    repop.s1, repop.s2, depop.r, abx.r, abx.s, timestep)