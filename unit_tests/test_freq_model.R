# Regression testing for binary model
rm(list=ls()) # Clean working environment
source("model_frequency.R") # Load model for testing

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


############################################# Function tests ##################################################
#################################### Test patient matrix generation ####################################

# Cases: Single time step
test_mean <- 5
tolerance <- 1
patient_mat.s <- patient.table(n.bed = 20, n.day = 300, mean.max.los = test_mean, timestep=1)
# Expected output: At high points, patient table should give exponential distribution with mean given within tolerance
hist(table(patient_mat.s)) # eyeball that this looks like an exponential distribution and that maximum value at tail makes sense
stopifnot(abs(mean(table(patient_mat.s)) - test_mean) < tolerance) # Make sure gives correct mean
stopifnot(dim(patient_mat.s) == c(300, 20)) # Make sure dimensions are correct

# Cases: Multiple time step
tolerance <- 1*3
patient_mat.m <- patient.table(n.bed = 100, n.day = 90, mean.max.los = test_mean, timestep=3)
# Expected output: 
hist(table(patient_mat.m)) # eyeball that this looks like an exponential distribution and that maximum value at tail makes sense
stopifnot(abs(mean(table(patient_mat.m)) - test_mean*3) < tolerance) # Make sure gives correct mean
stopifnot(dim(patient_mat.m) == c(270, 100)) # Make sure dimensions are correct

#################################### Test los.array summary ####################################

# Cases: Single time step
los_duration.s <- summary.los(patient_mat.s)
# Expected output: 
hist(los_duration.s[2, ]) # check distribution of length of stay, should be exponential
hist(table(patient_mat.s)) # chould be the same as input, los.array should = table(patient.table.output)
stopifnot(max(patient_mat.s) == max(los_duration.s)) # highest number from array should be exactly the same as highest id in patient.matrix

# Cases: Multiple timesteps
los_duration.m <- summary.los(patient_mat.m)
# Expected output: 
hist(los_duration.m[2, ]) # check distribution of length of stay, should be exponential
hist(table(patient_mat.m)) # chould be the same as input, los.array should = table(patient.table.output)
stopifnot(max(patient_mat.m) == max(los_duration.m)) # highest number from array should be exactly the same as highest id in patient.matrix

#################################### Test abx matrix generation ####################################

# Cases: Single time step
test_p <- 0.3
test_mean <- 10
tolerance <- 1
# Should be equivalent to simple_model test
abx.matrix.s <- abx.table(patient_mat.s, los_duration.s, p.s=test_p, p.r.day1=0, p.r.dayafter=0,
                          meanDur.s=test_p, meanDur.r=0, sdDur=1, timestep=1)
# Expected output: 
# To use commented test below, need to edit code: c(rep(1, abx_person), rep(0, max_days-abx_person)) to c(rep(1, abx_person), rep(1, max_days-abx_person))
# stopifnot(abs(sum(abx.matrix.s > 0)/length(abx.matrix.s) - test_p) < tolerance) # overall number of 1s approx. probability

# 100% effectiveness, no clipping of antibiotics
test_p <- 1
test_mean <- 10
tolerance <- 1
patient_mat.s <- patient.table(n.bed = 10, n.day = 3000, mean.max.los = 100, timestep=1) # have days long enough to not create clipping
los_duration.s <- summary.los(patient_mat.s)
abx.matrix.s <- abx.table(patient_mat.s, los_duration.s, p.s=test_p, p.r.day1=0, p.r.dayafter=0,
                          meanDur.s=test_mean, meanDur.r=0, sdDur=1, timestep=1)
abx_summary <- table(abx.matrix.s*patient_mat.s)
# Expected output: 
hist(abx_summary[2:length(abx_summary)], breaks=20) # distribution of abx duration shape should be truncated norm as used
stopifnot(abs(mean(abx_summary[2:length(abx_summary)]) - test_mean) < tolerance) # mean of abx duration is within tolerance
stopifnot(dim(abx.matrix.s) == dim(patient_mat.s)) # dimensions must equal patient.matrix

# 100% effectiveness, all clipping of antibiotics
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

# Cases: Multiple time step
test_p <- 0.3
test_mean <- 10
tolerance <- 1
abx.matrix.s <- abx.table(patient_mat.s, los_duration.s, p.s=test_p, p.r.day1=0, p.r.dayafter=0,
                          meanDur.s=test_mean, meanDur.r=0, sdDur=1, timestep=2)
# Expected output: 
# To use commented test below, need to edit code: c(rep(1, abx_person), rep(0, max_days-abx_person)) to c(rep(1, abx_person), rep(1, max_days-abx_person))
# stopifnot(abs(sum(abx.matrix.s > 0)/length(abx.matrix.s) - test_p) < tolerance) # overall number of 1s approx. probability

# 100% effectiveness, no clipping of antibiotics
test_p <- 1
test_mean <- 10
tolerance <- 1*2
patient_mat.s <- patient.table(n.bed = 20, n.day = 3000, mean.max.los = 100, timestep=2) # have days long enough to not create clipping
los_duration.s <- summary.los(patient_mat.s)
abx.matrix.s <- abx.table(patient_mat.s, los_duration.s, p.s=test_p, p.r.day1=0, p.r.dayafter=0,
                          meanDur.s=test_mean, meanDur.r=0, sdDur=1, timestep=2)
abx_summary <- table(abx.matrix.s*patient_mat.s)
# Expected output: 
hist(abx_summary[2:length(abx_summary)], breaks=20) # distribution of abx duration shape should be truncated norm as used
stopifnot(abs(mean(abx_summary[2:length(abx_summary)]) - test_mean*2) < tolerance) # mean of abx duration is within tolerance
stopifnot(dim(abx.matrix.s) == dim(patient_mat.s)) # dimensions must equal patient.matrix

# 100% effectiveness, all clipping of antibiotics
# should track exponential function of length of stay
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

# Case: Single timestep, multiple abx prob
# Case: Multiple timestep, multiple abx prob
# Case: Single timestep, check prob of starting abx per day
# Case: Multiple timestep, check prob of starting abx per day

#################################### Test starting bacteria generation ####################################

# Case: single values
# Default values
t_mean <- 4.0826
t_sd <- 1.1218
r_mean <- 1.7031
r_sd <- 1.8921
colo.matrix <- colo.table(patient_mat.s, los_duration.s, t_mean, t_sd, r_mean, r_sd)
stopifnot(sum(is.na(colo.matrix[[1]][1, ])) == 0) # Top row is fully initialized, S
stopifnot(sum(is.na(colo.matrix[[2]][1, ])) == 0) # Top row is fully initialized, R
stopifnot(dim(colo.matrix[[1]]) == dim(patient_mat.s)) # dimension must equal patient.matrix, S
stopifnot(dim(colo.matrix[[2]]) == dim(patient_mat.s)) # dimension must equal patient.matrix, R
patient_idx <- cumsum(c(1, los_duration.s[2,]))
stopifnot(sum(!(which(!is.na(colo.matrix[[1]])) == patient_idx[-length(patient_idx)])) == 0) # starts on same position as patient id, S
patient_idx <- cumsum(c(1, los_duration.s[2,]))
stopifnot(sum(!(which(!is.na(colo.matrix[[2]])) == patient_idx[-length(patient_idx)])) == 0) # starts on same position as patient id, R

############################################ Integration tests ##################################################
# Test diff_prevalence
diff_prevalence(n.bed=20, mean.max.los=5, p.s=0.10, p.r.day1=0.10, p.r.dayafter=0.10,
                K=1000, t_mean=4.0826, t_sd=1.1218, r_mean=1.7031, r_sd=1.8921,
                pi_r=0.1, r_thres=10, r_growth=2, r_trans=10, 
                abxr_killr=5, abxr_kills=5, abxs_kills=5,
                short_dur.s=4, long_dur.s=14, short_dur.r = 2, long_dur.r = 5, sdDur = 5)

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