# Regression testing for binary model
rm(list=ls()) # Clean working environment
source("model_binary.R") # Load model for testing

# Please add as appropriate, be as pendantic as much as possible
# Cleanest way to run. Source above, select only block of test case, try running

############################################# Function tests ##################################################
########### Test patient matrix generation

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


########### Test los.array summary

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


########### Test abx matrix generation

# Cases: Single time step
test_p <- 0.3
test_mean <- 10
tolerance <- 1
abx.matrix.s <- abx.table(patient_mat.s, los_duration.s, p=test_p, meanDur=test_mean, sdDur=1, timestep = 1)
# Expected output: 
# To use commented test below, need to edit code: c(rep(1, abx_person), rep(0, max_days-abx_person)) to c(rep(1, abx_person), rep(1, max_days-abx_person))
# stopifnot(abs(sum(abx.matrix.s > 0)/length(abx.matrix.s) - test_p) < tolerance) # overall number of 1s approx. probability

# 100% effectiveness, no clipping of antibiotics
test_p <- 1
test_mean <- 10
tolerance <- 1
patient_mat.s <- patient.table(n.bed = 20, n.day = 3000, mean.max.los = 100, timestep=1) # have days long enough to not create clipping
los_duration.s <- summary.los(patient_mat.s)
abx.matrix.s <- abx.table(patient_mat.s, los_duration.s, p=test_p, meanDur=test_mean, sdDur=1, timestep=1)
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
patient_mat.s <- patient.table(n.bed = 20, n.day = 30000, mean.max.los = test_los_mean, timestep=1) # have days short to create clipping
los_duration.s <- summary.los(patient_mat.s)
abx.matrix.s <- abx.table(patient_mat.s, los_duration.s, p=test_p, meanDur=test_mean, sdDur=1, timestep=1)
abx_summary <- table(abx.matrix.s*patient_mat.s)
# Expected output: 
hist(abx_summary[2:length(abx_summary)], breaks=20) # distribution of abx duration shape should be exponential function, small bump at test_mean for "lucky hits"
stopifnot(abs(mean(abx_summary[2:length(abx_summary)]) - test_los_mean) < tolerance) # mean of abx duration is mean of los instead
stopifnot(dim(abx.matrix.s) == dim(patient_mat.s)) # dimensions must equal patient.matrix

# Cases: Multiple time step
test_p <- 0.3
test_mean <- 10
tolerance <- 1
abx.matrix.s <- abx.table(patient_mat.s, los_duration.s, p=test_p, meanDur=test_mean, sdDur=1, timestep = 1)
# Expected output: 
# To use commented test below, need to edit code: c(rep(1, abx_person), rep(0, max_days-abx_person)) to c(rep(1, abx_person), rep(1, max_days-abx_person))
# stopifnot(abs(sum(abx.matrix.s > 0)/length(abx.matrix.s) - test_p) < tolerance) # overall number of 1s approx. probability

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

