# Regression testing for simple model
rm(list=ls()) # Clean working environment
source("model_simple.R") # Load model for testing

# Please add as appropriate, be as pendantic as much as possible
# Cleanest way to run. Source above lines, select only block of test case, try running

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


############################################# Function tests ##################################################
#################################### Test patient matrix generation ####################################
# ---------------------------------------------------------------------------------------
# Test name: simple run
# Test summary: Tests patient.table functions exits with no errors, gives correct mean, and gives correct shape
# Function tested: patient.table
# Pre-conditions:
rm(list=ls()) # Clean working environment
source("model_simple.R") # Load model for testing
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
source("model_simple.R") # Load model for testing
test_mean <- 7
tolerance <- 1*3
patient_mat.m <- patient.table(n.bed = 100, n.day = 90, mean.max.los = test_mean, timestep=3)
# Expected result: looks like an exponential distribution and that maximum value at tail makes sense
hist(table(patient_mat.m))
# Expected result: mean from table is the same mean given within tolerance
stopifnot(abs(mean(table(patient_mat.m)) - test_mean*3) < tolerance)
# Expected result: dimensions are correct dimensions
stopifnot(dim(patient_mat.m) == c(270, 100))

#################################### Test los.array summary ####################################
# ---------------------------------------------------------------------------------------
# Test name: 
# Test summary: 
# Function tested: summary.los
# Pre-conditions: 
rm(list=ls()) # Clean working environment
source("model_simple.R") # Load model for testing
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
source("model_simple.R") # Load model for testing
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
# Test name: Single time step basic run
# Test summary:
# Function tested: abx.table
# Pre-conditions: 
rm(list=ls()) # Clean working environment
source("model_simple.R") # Load model for testing
test_p <- 0.3
test_mean <- 10
tolerance <- 1
patient_mat.s <- patient.table(n.bed = 100, n.day = 90, mean.max.los = test_mean, timestep=1)
los_duration.s <- summary.los(patient_mat.s)
# Expected result: 
abx.matrix.s <- abx.table(patient_mat.s, los_duration.s, p=test_p, meanDur=test_mean, sdDur=1, timestep = 1)
# Expected output: 
# To use commented test below, need to edit code: c(rep(1, abx_person), rep(0, max_days-abx_person)) to c(rep(1, abx_person), rep(1, max_days-abx_person))
# stopifnot(abs(sum(abx.matrix.s > 0)/length(abx.matrix.s) - test_p) < tolerance) # overall number of 1s approx. probability

# ---------------------------------------------------------------------------------------
# Test name: p, meanDur
# Test summary: # 100% effectiveness, no clipping of antibiotics
# Function tested: abx.table
# Pre-conditions: 
rm(list=ls()) # Clean working environment
source("model_simple.R") # Load model for testing
test_p <- 1
test_mean <- 10
tolerance <- 1
patient_mat.s <- patient.table(n.bed = 20, n.day = 3000, mean.max.los = 100, timestep=1) # have days long enough to not create clipping
los_duration.s <- summary.los(patient_mat.s)
abx.matrix.s <- abx.table(patient_mat.s, los_duration.s, p=test_p, meanDur=test_mean, sdDur=1, timestep=1)
abx_summary <- table(abx.matrix.s*patient_mat.s)
# Expected result: # distribution of abx duration shape should be truncated norm as used
hist(abx_summary[2:length(abx_summary)], breaks=20) 
# Expected result: # mean of abx duration is within tolerance
stopifnot(abs(mean(abx_summary[2:length(abx_summary)]) - test_mean) < tolerance)
# Expected result: # dimensions must equal patient.matrix
stopifnot(dim(abx.matrix.s) == dim(patient_mat.s)) 

# ---------------------------------------------------------------------------------------
# Test name: meanDur > mean.max.los
# Test summary: # 100% effectiveness, all clipping of antibiotics
# Function tested: abx.table
# Pre-conditions: 
test_p <- 1
test_mean <- 10
test_los_mean <- 5 # have days short to create clipping
tolerance <- 1
patient_mat.s <- patient.table(n.bed = 20, n.day = 30000, mean.max.los = test_los_mean, timestep=1) 
los_duration.s <- summary.los(patient_mat.s)
abx.matrix.s <- abx.table(patient_mat.s, los_duration.s, p=test_p, meanDur=test_mean, sdDur=1, timestep=1)
abx_summary <- table(abx.matrix.s*patient_mat.s)
# Expected result: # should track exponential function of length of stay
# distribution of abx duration shape should be exponential function, small bump at test_mean for "lucky hits"
hist(abx_summary[2:length(abx_summary)], breaks=20) 
# Expected result: # mean of abx duration is mean of los instead
stopifnot(abs(mean(abx_summary[2:length(abx_summary)]) - test_los_mean) < tolerance) 
# Expected result: # dimensions must equal patient.matrix
stopifnot(dim(abx.matrix.s) == dim(patient_mat.s)) 

# ---------------------------------------------------------------------------------------
# Test name: 
# Test summary: 
# Function tested: 
# Pre-conditions: 
# Expected result: 
# Cases: Multiple time step
test_p <- 0.3
test_mean <- 10
tolerance <- 1
abx.matrix.s <- abx.table(patient_mat.s, los_duration.s, p=test_p, meanDur=test_mean, sdDur=1, timestep = 1)
# Expected output: 
# To use commented test below, need to edit code: c(rep(1, abx_person), rep(0, max_days-abx_person)) to c(rep(1, abx_person), rep(1, max_days-abx_person))
# stopifnot(abs(sum(abx.matrix.s > 0)/length(abx.matrix.s) - test_p) < tolerance) # overall number of 1s approx. probability

# ---------------------------------------------------------------------------------------
# Test name: 
# Test summary: 
# Function tested: 
# Pre-conditions: 
# Expected result: 
# 100% effectiveness, no clipping of antibiotics
test_p <- 1
test_mean <- 10
tolerance <- 1*2
patient_mat.s <- patient.table(n.bed = 20, n.day = 3000, mean.max.los = 100, timestep=2) # have days long enough to not create clipping
los_duration.s <- summary.los(patient_mat.s)
abx.matrix.s <- abx.table(patient_mat.s, los_duration.s, p=test_p, meanDur=test_mean, sdDur=1, timestep=2)
abx_summary <- table(abx.matrix.s*patient_mat.s)
# Expected output: 
hist(abx_summary[2:length(abx_summary)], breaks=20) # distribution of abx duration shape should be truncated norm as used
stopifnot(abs(mean(abx_summary[2:length(abx_summary)]) - test_mean*2) < tolerance) # mean of abx duration is within tolerance
stopifnot(dim(abx.matrix.s) == dim(patient_mat.s)) # dimensions must equal patient.matrix

# ---------------------------------------------------------------------------------------
# Test name: 
# Test summary: 
# Function tested: 
# Pre-conditions: 
# Expected result: 
# 100% effectiveness, all clipping of antibiotics
# should track exponential function of length of stay
test_p <- 1
test_mean <- 10
test_los_mean <- 5 
tolerance <- 1*2
patient_mat.s <- patient.table(n.bed = 20, n.day = 3000, mean.max.los = test_los_mean, timestep=2) # have days short to create clipping
los_duration.s <- summary.los(patient_mat.s)
abx.matrix.s <- abx.table(patient_mat.s, los_duration.s, p=test_p, meanDur=test_mean, sdDur=1, timestep=2)
abx_summary <- table(abx.matrix.s*patient_mat.s)
# Expected output: 
hist(abx_summary[2:length(abx_summary)], breaks=20) # distribution of abx duration shape should be exponential function, small bump at test_mean for "lucky hits"
stopifnot(abs(mean(abx_summary[2:length(abx_summary)]) - test_los_mean*2) < tolerance) # mean of abx duration is mean of los instead
stopifnot(dim(abx.matrix.s) == dim(patient_mat.s)) # dimensions must equal patient.matrix


#################################### Test starting bacteria generation ####################################
# ---------------------------------------------------------------------------------------
# Test name: 
# Test summary: 
# Function tested: 
# Pre-conditions: 
# Expected result: 
# Cases: Single time step
prob_StartBact_R <- 0.2
prop_S_nonR <- 0.3
tolerance <- 0.1
patient.matrix <- patient.table(n.bed = 20, n.day = 300, mean.max.los = 5, timestep=1)
los.array <- summary.los(patient.matrix)
colo.matrix <- colo.table(patient.matrix, los.array, prob_StartBact_R, prop_S_nonR)
# Expected output: 
stopifnot(abs(sum(colo.matrix == "R", na.rm=TRUE)/sum(!is.na(colo.matrix)) - prob_StartBact_R) < tolerance) # generated R probs is correct
stopifnot(abs(sum(colo.matrix == "S", na.rm=TRUE)/sum(!is.na(colo.matrix)) - prop_S_nonR*(1-prob_StartBact_R)) < tolerance) # generated R probs is correct
stopifnot(abs(sum(colo.matrix == "ss", na.rm=TRUE)/sum(!is.na(colo.matrix)) - ((1-prop_S_nonR)*(1-prob_StartBact_R))) < tolerance) # generated R probs is correct
stopifnot(sum(is.na(colo.matrix[1, ])) == 0) # Top row is fully initialized
stopifnot(dim(colo.matrix) == dim(patient.matrix)) # dimension must equal patient.matrix
patient_idx <- cumsum(c(1, los.array[2,]))
stopifnot(sum(!(which(!is.na(colo.matrix)) == patient_idx[-length(patient_idx)])) == 0) # starts on same position as patient id

# ---------------------------------------------------------------------------------------
# Test name: 
# Test summary: 
# Function tested: 
# Pre-conditions: 
# Expected result: 
# Cases: Multiple time step
prob_StartBact_R <- 0.6
prop_S_nonR <- 0.2
tolerance <- 0.1
patient.matrix <- patient.table(n.bed = 20, n.day = 300, mean.max.los = 5, timestep=3)
los.array <- summary.los(patient.matrix)
colo.matrix <- colo.table(patient.matrix, los.array, prob_StartBact_R, prop_S_nonR)
# Expected output: 
stopifnot(abs(sum(colo.matrix == "R", na.rm=TRUE)/sum(!is.na(colo.matrix)) - prob_StartBact_R) < tolerance) # generated R probs is correct
stopifnot(abs(sum(colo.matrix == "S", na.rm=TRUE)/sum(!is.na(colo.matrix)) - prop_S_nonR*(1-prob_StartBact_R)) < tolerance) # generated R probs is correct
stopifnot(abs(sum(colo.matrix == "ss", na.rm=TRUE)/sum(!is.na(colo.matrix)) - ((1-prop_S_nonR)*(1-prob_StartBact_R))) < tolerance) # generated R probs is correct
stopifnot(sum(is.na(colo.matrix[1, ])) == 0) # Top row is fully initialized
stopifnot(dim(colo.matrix) == dim(patient.matrix))# dimension must equal patient.matrix
patient_idx <- cumsum(c(1, los.array[2,]))
stopifnot(sum(!(which(!is.na(colo.matrix)) == patient_idx[-length(patient_idx)])) == 0) # starts on same position as patient id


#################################### Test daily updates ####################################

# ---------------------------------------------------------------------------------------
# Test name: 
# Test summary: 
# Function tested: 
# Pre-conditions: 
# Expected result: 
# Cases: Single time step
tolerance <- 0.02
pi_ssr <- 0.05
patient.matrix <- patient.table(n.bed = 20, n.day = 300, mean.max.los = 3, timestep=1)
los.array <- summary.los(patient.matrix)
abx.matrix <- abx.table(patient.matrix, los.array, p=0.3, meanDur=5, sdDur=1, timestep=1)
colo.matrix <- colo.table(patient.matrix, los.array, prob_StartBact_R=0.1, prop_S_nonR=0.5)
update <- nextDay(patient.matrix, los.array, abx.matrix, colo.matrix, 
                    bif=1, pi_ssr, repop.s1=0, mu_r=0, abx.clear=1)
# Expected output:
stopifnot(sum(is.na(update)) == 0) # all days were filled
stopifnot(sum(!(update == "R" | update == "S" | update == "ss")) == 0) # all days had expected values
colo_idx <- which(!is.na(colo.matrix))
stopifnot(sum(!colo.matrix[colo_idx] == update[colo_idx]) == 0) # starting conditions on colo.matrix was not overwritten

# ---------------------------------------------------------------------------------------
# Test name: 
# Test summary: 
# Function tested: 
# Pre-conditions: 
# Expected result: 
# update probability R -> S holds (mu_r)
# No abx, no transmission, 
tolerance <- 0.02
mu_r <- 0.7
patient.matrix <- patient.table(n.bed = 20, n.day = 3000, mean.max.los=3, timestep=1)
los.array <- summary.los(patient.matrix)
abx.matrix <- abx.table(patient.matrix, los.array, p=0, meanDur=5, sdDur=1, timestep=1)
colo.matrix <- colo.table(patient.matrix, los.array, prob_StartBact_R=0.6, prop_S_nonR=0)
update <- nextDay(patient.matrix, los.array, abx.matrix, colo.matrix,
                  bif=0, pi_ssr=0, repop.s1=0, mu_r, abx.clear=1)
num_update <- matrix(NA, nrow=nrow(update), ncol=ncol(update))
num_update[update == "ss"] <- 0
num_update[update == "R"] <- 5
num_update[update == "S"] <- 6
parse_list <- split(num_update, patient.matrix) 
state_change <- sum(unlist(lapply(parse_list, function(x) diff(x))) == 1)
count_update <- matrix(NA, nrow=nrow(update), ncol=ncol(update))
count_update[update == "ss"] <- 0
count_update[update == "R"] <- 1
count_update[update == "S"] <- 0
parse_list <- split(count_update, patient.matrix)
r_count <- sum(unlist(lapply(parse_list, function(x) x[-length(x)]))) # count only R not at the end of each patient
# Expect output:
stopifnot(abs(state_change/r_count - mu_r) < tolerance) # update probability R -> S holds (mu_r)

# ---------------------------------------------------------------------------------------
# Test name: 
# Test summary: 
# Function tested: 
# Pre-conditions: 
# Expected result: 
# update probability ss -> S holds (repop.s1)
tolerance <- 0.02
repop.s1 <- 0.5
patient.matrix <- patient.table(n.bed = 20, n.day = 3000, mean.max.los=3, timestep=1)
los.array <- summary.los(patient.matrix)
abx.matrix <- abx.table(patient.matrix, los.array, p=0, meanDur=5, sdDur=1, timestep=1)
colo.matrix <- colo.table(patient.matrix, los.array, prob_StartBact_R=0.2, prop_S_nonR=0.5)
update <- nextDay(patient.matrix, los.array, abx.matrix, colo.matrix,
                  bif=0, pi_ssr=0, repop.s1, mu_r=0, abx.clear=1)
num_update <- matrix(NA, nrow=nrow(update), ncol=ncol(update))
num_update[update == "ss"] <- 5
num_update[update == "R"] <- 0
num_update[update == "S"] <- 6
parse_list <- split(num_update, patient.matrix) 
state_change <- sum(unlist(lapply(parse_list, function(x) diff(x))) == 1)
count_update <- matrix(NA, nrow=nrow(update), ncol=ncol(update))
count_update[update == "ss"] <- 1
count_update[update == "R"] <- 0
count_update[update == "S"] <- 0
parse_list <- split(count_update, patient.matrix)
ss_count <- sum(unlist(lapply(parse_list, function(x) x[-length(x)]))) # count only ss not at the end of each patient
# Expect output:
stopifnot(abs(state_change/ss_count - repop.s1) < tolerance) # update probability ss -> S holds (repop.s1)

# ---------------------------------------------------------------------------------------
# Test name: 
# Test summary: 
# Function tested: 
# Pre-conditions: 
# Expected result: 
# update probability S -> ss holds (p*abx.clear)
# abx.clear = 1, all patients have the same length of stay, each antibiotic given based on flat probability
# remove effect of patient exponential distribution and abx truncated norm distribution, to check only update works
tolerance <- 0.02
p <- 0.3
# generate patient matrix where each person is equally given 3 days
patient.matrix <- patient.table.equal(20, 300, 3)
los.array <- summary.los(patient.matrix)
roll <- runif(100*30*1, 0, 1)
abx.matrix <- matrix(as.numeric(roll < p), nrow=nrow(patient.matrix), ncol=ncol(patient.matrix))
# sum(abx.matrix == 1)/length(abx.matrix) # Check abx gives correct probability
colo.matrix <- colo.table(patient.matrix, los.array, prob_StartBact_R=0, prop_S_nonR=1)
update <- nextDay(patient.matrix, los.array, abx.matrix, colo.matrix,
                  bif=0, pi_ssr=0, repop.s1=0, mu_r=0, abx.clear=1)
num_update <- matrix(NA, nrow=nrow(update), ncol=ncol(update))
num_update[update == "ss"] <- 6
num_update[update == "R"] <- 0
num_update[update == "S"] <- 5
parse_list <- split(num_update, patient.matrix) 
state_change <- sum(unlist(lapply(parse_list, function(x) diff(x))) == 1)
count_update <- matrix(NA, nrow=nrow(update), ncol=ncol(update))
count_update[update == "ss"] <- 0
count_update[update == "R"] <- 0
count_update[update == "S"] <- 1
parse_list <- split(count_update, patient.matrix)
ss_count <- sum(unlist(lapply(parse_list, function(x) x[-length(x)])))  # count only S not at the end of each patient
# Expect output:
stopifnot(abs(state_change/ss_count - p) < tolerance) # update probability S -> ss holds (p*abx.clear)

# ---------------------------------------------------------------------------------------
# Test name: 
# Test summary: 
# Function tested: 
# Pre-conditions: 
# Expected result: 
# abx.clear = 0.36, check abx.clear augments p of getting abx correctly
tolerance <- 0.02
p <- 0.3
abx.clear <- 0.36
# generate patient matrix where each person is equally given 3 days
patient.matrix <- patient.table.equal(20, 300, 3)
los.array <- summary.los(patient.matrix)
roll <- runif(100*30*1, 0, 1)
abx.matrix <- matrix(as.numeric(roll < p), nrow=nrow(patient.matrix), ncol=ncol(patient.matrix))
# sum(abx.matrix == 1)/length(abx.matrix) # Check abx gives correct probability
colo.matrix <- colo.table(patient.matrix, los.array, prob_StartBact_R=0, prop_S_nonR=1)
update <- nextDay(patient.matrix, los.array, abx.matrix, colo.matrix,
                  bif=0, pi_ssr=0, repop.s1=0, mu_r=0, abx.clear)
num_update <- matrix(NA, nrow=nrow(update), ncol=ncol(update))
num_update[update == "ss"] <- 6
num_update[update == "R"] <- 0
num_update[update == "S"] <- 5
parse_list <- split(num_update, patient.matrix) 
state_change <- sum(unlist(lapply(parse_list, function(x) diff(x))) == 1)
count_update <- matrix(NA, nrow=nrow(update), ncol=ncol(update))
count_update[update == "ss"] <- 0
count_update[update == "R"] <- 0
count_update[update == "S"] <- 1
parse_list <- split(count_update, patient.matrix)
ss_count <- sum(unlist(lapply(parse_list, function(x) x[-length(x)])))  # count only S not at the end of each patient
# Expect output:
stopifnot(abs(state_change/ss_count - p*abx.clear) < tolerance) 

# ---------------------------------------------------------------------------------------
# Test name: pi_ssr
# Test summary: update probability ss -> R holds (pi_ssr)
# Function tested: nextDay 
# Pre-conditions: 
tolerance <- 0.02
patient.matrix<-patient.table(n.bed=20, n.day=50, mean.max.los=5, timestep=1)
los.array<- summary.los(patient.matrix=patient.matrix)
abx.matrix<- abx.table(patient.matrix=patient.matrix, los.array=los.array, p=0.2, meanDur=3, sdDur=1, timestep=1)
colo.matrix<- colo.table(patient.matrix=patient.matrix, los.array=los.array, prob_StartBact_R=0.5, prop_S_nonR=0.3)
# Expected result: 
bif=1
repop.s1=0
mu_r=0
abx.clear=0.5 
timestep=1
pi_ssr <- 0.3
pi_Sr <- pi_ssr - (bif*pi_ssr)

#Check S (pi_Sr)
for(i in 2:nrow(patient.matrix)){
# Get the previous row subset of the entire matrix which represents previous day or timestep
prev_step <- colo.matrix[i-1, ]
# Get the indices which are already filled by a patient entering the ward for exclusion
already_filled <- which(!is.na(colo.matrix[i, ]))
# Get the column indices which contain S in the previous day
S <- which(prev_step == "S")
# Remove column indices that already have a starting bacterial state filled in
S <- S[!(S %in% already_filled)]
# count if there are any S on the previous day
s_num <- length(S)
# if there is any S (number of S > 0) in the previous day
r_num<-5
if(s_num){
    # roll for transmission of R
    prob_r <- 1-((1-pi_Sr)^r_num)
    # Roll a random number for each R on the previous day for clearance
    roll_clear <- runif(s_num, 0, 1)
    # All column indices which roll < probability of clearance AND there is antibiotic used on that patient-timestep
    clear_idx <- S[abx.matrix[i, S] & (roll_clear < abx.clear)]
    # Clear those that pass roll and use abx to ss
    colo.matrix[i, clear_idx] <- "ss"
    # Removed those that have been cleared by abx from list of S indices
    S <- S[!(S %in% clear_idx)]
    # Roll a random number for each remaining S for chance of selection to
    roll_trans <- runif(length(S), 0, 1)
    r_idx <- S[roll_trans < prob_r]
    same_idx <- S[roll_trans >= prob_r]
    colo.matrix[i, r_idx] <- "R"
    colo.matrix[i, same_idx] <- "S"
}
}

update <- colo.matrix
num_update <- matrix(NA, nrow=nrow(update), ncol=ncol(update))
num_update[update == "ss"] <- 5
num_update[update == "R"] <- 6
num_update[update == "S"] <- 0
parse_list <- split(num_update, patient.matrix) 
state_change_R <- sum(unlist(lapply(parse_list, function(x) diff(x))) == 6, na.rm = T)
state_no_change <- sum(unlist(lapply(parse_list, function(x) diff(x))) == 0, na.rm = T)
state_change_ss <- sum(unlist(lapply(parse_list, function(x) diff(x))) == 5, na.rm = T)
state_change_ss/sum(state_no_change,state_change_ss, state_change_R)
prob.r<-state_change_R/sum(state_no_change,state_change_ss, state_change_R)
1-(1-prob.r)^1/5


count_update[update[, -(1:(r_prop*20))] == "ss"] <- 1
count_update[update[, -(1:(r_prop*20))] == "R"] <- 0
count_update[update[, -(1:(r_prop*20))] == "S"] <- 0
parse_list <- split(count_update, patient.matrix[, -(1:(r_prop*20))])
ss_count <- sum(unlist(lapply(parse_list, function(x) x[-length(x)])))
# Expected output:
state_change/ss_count
stopifnot(abs(state_change/length(colo.matrix == "ss") - pi_ssr) < tolerance) # update probability ss -> R holds (pi_ssr)

# ---------------------------------------------------------------------------------------
# Test name: 
# Test summary: 
# Function tested: 
# Pre-conditions: 
# Expected result: 
# No abx, no starting R, get probability of transmission from S to R, (pi_Sr)
#stopifnot(abs(state_change/(length(colo.matrix) - length(colo_idx)) - pi_ssr) < tolerance) # update probability S -> R holds (pi_Sr)


# ---------------------------------------------------------------------------------------
# Test name: 
# Test summary: 
# Function tested: 
# Pre-conditions: 
# Expected result: 
# Cases: Multiple time steps

############################################ Integration tests ##################################################
# Test diff_prevalence
diff_prevalence(n.bed=20, mean.max.los=3,
                prob_StartBact_R=0.2, prop_S_nonR=0.1,
                bif=0, pi_ssr=0.02, repop.s1=0.01, mu_r=0.01, abx.clear=0.3,
                p=0.4, short_dur=4, long_dur=10, sdDur = 5)

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
