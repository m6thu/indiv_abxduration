# Regression testing for simple model
source("model_simple.R")

# Please add as appropriate, be as pendantic as much as possible
# Cleanest way to run. Source above, select only block of test case, try running

############################################# Function tests ##################################################
########### Test patient matrix generation

# Cases: Single time step
test_mean <- 5
tolerance <- 1
patient_mat.s <- patient.table(n.bed = 20, n.day = 10, mean.max.los = test_mean, timestep=1)
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
tolerance <- 1
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
tolerance <- 1
patient_mat.s <- patient.table(n.bed = 20, n.day = 3000, mean.max.los = test_los_mean, timestep=2) # have days short to create clipping
los_duration.s <- summary.los(patient_mat.s)
abx.matrix.s <- abx.table(patient_mat.s, los_duration.s, p=test_p, meanDur=test_mean, sdDur=1, timestep=2)
abx_summary <- table(abx.matrix.s*patient_mat.s)
# Expected output: 
hist(abx_summary[2:length(abx_summary)], breaks=20) # distribution of abx duration shape should be exponential function, small bump at test_mean for "lucky hits"
stopifnot(abs(mean(abx_summary[2:length(abx_summary)]) - test_los_mean*2) < tolerance) # mean of abx duration is mean of los instead
stopifnot(dim(abx.matrix.s) == dim(patient_mat.s)) # dimensions must equal patient.matrix


########### Test starting bacteria generation

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


########### Test daily updates

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
num_update <- matrix(NA, nrow=nrow(update), ncol=ncol(update))
num_update[update == "ss"] <- 0
num_update[update == "R"] <- 1
num_update[update == "S"] <- 0
parse_list <- split(num_update, patient.matrix)
r_count <- sum(unlist(lapply(parse_list, function(x) x[-length(x)]))) # count only R not at the end of each patient
# Expect output:
stopifnot(abs(state_change/r_count - mu_r) < tolerance) # update probability R -> S holds (mu_r)

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
num_update <- matrix(NA, nrow=nrow(update), ncol=ncol(update))
num_update[update == "ss"] <- 1
num_update[update == "R"] <- 0
num_update[update == "S"] <- 0
parse_list <- split(num_update, patient.matrix)
ss_count <- sum(unlist(lapply(parse_list, function(x) x[-length(x)]))) # count only ss not at the end of each patient
# Expect output:
state_change/ss_count
stopifnot(abs(state_change/ss_count - repop.s1) < tolerance) # update probability ss -> S holds (repop.s1)

# update probability S -> ss holds (p*abx.clear)
# abx.clear = 1
tolerance <- 0.02
p <- 0.2
patient.matrix <- patient.table(n.bed = 20, n.day = 3000, mean.max.los=3, timestep=1)
los.array <- summary.los(patient.matrix)
abx.matrix <- abx.table(patient.matrix, los.array, p=p, meanDur=5, sdDur=1, timestep=1)
colo.matrix <- colo.table(patient.matrix, los.array, prob_StartBact_R=0.2, prop_S_nonR=0.5)
update <- nextDay(patient.matrix, los.array, abx.matrix, colo.matrix,
                  bif=0, pi_ssr=0, repop.s1, mu_r=0, abx.clear=1)
num_update <- matrix(NA, nrow=nrow(update), ncol=ncol(update))
num_update[update == "ss"] <- 5
num_update[update == "R"] <- 0
num_update[update == "S"] <- 6
parse_list <- split(num_update, patient.matrix) 
state_change <- sum(unlist(lapply(parse_list, function(x) diff(x))) == 1)
num_update <- matrix(NA, nrow=nrow(update), ncol=ncol(update))
num_update[update == "ss"] <- 1
num_update[update == "R"] <- 0
num_update[update == "S"] <- 0
parse_list <- split(num_update, patient.matrix)
ss_count <- sum(unlist(lapply(parse_list, function(x) x[-length(x)])))  # count only S not at the end of each patient
# Expect output:
state_change/ss_count
stopifnot(abs(state_change/ss_count - repop.s1) < tolerance) # update probability S -> ss holds (p*abx.clear)

# abx.clear = 0.6

# update probability ss -> R holds (pi_ssr)
tolerance <- 0.02
pi_ssr <- 0.1
patient.matrix <- patient.table(n.bed = 20, n.day = 300, mean.max.los = 3, timestep=1)
los.array <- summary.los(patient.matrix)
abx.matrix <- abx.table(patient.matrix, los.array, p=0.3, meanDur=5, sdDur=1, timestep=1)
colo.matrix <- colo.table(patient.matrix, los.array, prob_StartBact_R=0.1, prop_S_nonR=0.5)
update <- nextDay(patient.matrix, los.array, abx.matrix, colo.matrix, 
                  bif=1, pi_ssr, repop.s1=0, mu_r=0, abx.clear=1)
num_update <- matrix(NA, nrow=nrow(update), ncol=ncol(update))
num_update[update == "ss"] <- 5
num_update[update == "R"] <- 6
num_update[update == "S"] <- 1
parse_list <- split(num_update, patient.matrix) 
state_change <- sum(unlist(lapply(parse_list, function(x) diff(x))) == 1)
colo_idx <- which(colo.matrix == "ss")
state_change/(sum(update == "ss"))
stopifnot(abs(state_change/length(colo.matrix == "ss") - pi_ssr) < tolerance) # update probability ss -> R holds (pi_ssr)

# No abx, no starting R, get probability of transmission from S to R, (pi_Sr)
tolerance <- 0.02
pi_ssr <- 0.5 # bif = 1, therefore pi_Sr = pi_ssr
patient.matrix <- patient.table(n.bed = 20, n.day = 300, mean.max.los=3, timestep=1)
los.array <- summary.los(patient.matrix)
abx.matrix <- abx.table(patient.matrix, los.array, p=0, meanDur=5, sdDur=1, timestep=1)
colo.matrix <- colo.table(patient.matrix, los.array, prob_StartBact_R=0.1, prop_S_nonR=0.5)
update <- nextDay(patient.matrix, los.array, abx.matrix, colo.matrix,
                  bif=1, pi_ssr, repop.s1=0, mu_r=0, abx.clear=1)
num_update <- matrix(NA, nrow=nrow(update), ncol=ncol(update))
num_update[update == "ss"] <- 0
num_update[update == "S"] <- 5
num_update[update == "R"] <- 6
parse_list <- split(num_update, patient.matrix) 
state_change <- sum(unlist(lapply(parse_list, function(x) diff(x))) == 1)
colo_idx <- which(!is.na(colo.matrix))
# Expected output:
state_change/sum(colo.matrix == "S")
state_change/(length(colo.matrix) - length(colo_idx))
#stopifnot(abs(state_change/(length(colo.matrix) - length(colo_idx)) - pi_ssr) < tolerance) # update probability S -> R holds (pi_Sr)


############################################ Integration tests ##################################################
# Test diff_prevalence
diff_prevalence(n.bed=20, mean.max.los=3, timestep=1,
                prob_StartBact_R=0.4, prop_S_nonR=0.3,
                bif=0.003, pi_ssr=0.03, repop.s1=0, mu_r=0, abx.clear=0.1,
                p=0.1, short_dur=4, long_dur=10, sdDur = 5)

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
