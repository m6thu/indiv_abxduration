#######Modelling Day Project######
#######Parameter exploration######
################################### Dependencies and functions ################################################

# SAMPLE PARAMETER SPACE 
# load libraries 
require(pse) #load pse package for Latin Hypercube
require(sensitivity) #load sensitivity package for sensitivity analysis
require(parallel) # load parallel processing package to use multiple cores on computer (or cluster)
require(Rcpp) #optimising functions

cl <- makeCluster(detectCores())

model <- 'binary'
#source(paste0("model_binary.R"))

clusterCall(cl, function() {source('model_binary.R')})
# source functions on all cores
modelRun.binary <- function (data.df) { #data.df is a dataframe of the parameter values in columns 
    return(mapply(diff_prevalence, 
                  data.df[,1], data.df[,2], data.df[,3], data.df[,4], data.df[,5], 
                  data.df[,6], data.df[,7], data.df[,8], data.df[,9], 
                  data.df[,10], data.df[,11], data.df[,12], data.df[,13], data.df[,14], data.df[,15], 
                  data.df[,16], data.df[,17], data.df[,18],data.df[,19],
                  data.df[,20], data.df[,21], data.df[,22], data.df[,23], data.df[,24]
                  
    ))
}

################################## Define parameters and run LHS ####################################
#list parameters, the probability density functions from which the parameter values will be calculated, and what are the arguments to these functions
#list parameters together with name, so they are "linked" or not easily confused
parameters <- list(
    c("qunif", list(min=3, max=50), "n.bed"),              #n.bed; number of beds in the ward
    c("qunif", list(min=3, max=30), "mean.max.los"),       #mean.max.los; mean of length of stay (exponential distribution)
    c("qunif", list(min=0.1, max=0.5), "p.s"),             #probability of being prescribed narrow spectrum antibiotic
    c("qunif", list(min=0, max=0.4), "p.r.day1"),          #probability of being prescribed broad spectrum antibiotic on day 1 of admission 
    c("qunif", list(min=0, max= 0.2), "p.r.dayafter"),     #probability of being prescribed broad spectrum antibiotic after admission 
    c("qunif", list(min=0, max=0.9), "prob_StartBact_R"),  #probability of initial carriage of resistant organisms
    c("qunif", list(min=0, max=1), "prop_S_nonR"),         #proportion of S in (S+s): prob_start_S <- prop_S_nonR*(1-prob_StartBact_R)
    c("qunif", list(min=0, max=1), "prop_Sr_inR"),         #proportion of Sr in (r+R): prob_start_Sr <- prop_Sr_inR*prob_StartBact_R
    c("qunif", list(min=0, max=1), "prop_sr_inR"),         #proportion of sr in (r+r): prob_start_sr <- prop_sr_inR*prob_StartBact_R
    c("qunif", list(min=0, max=0.05), "pi_r2"),            #probability of being transmitted r to ss (ss—> ssr)
    c("qunif", list(min=0, max=1), "bif"),                 #bacterial interference factor (pi_r2 = pi_r1 * bif )
    c("qunif", list(min=0, max=0.005), "mu1"),             #probability of being decolonised to S (Sr—> S) 
    c("qunif", list(min=0, max=0.005), "mu2"),             #probability of being decolonised to S (sr—> s) 
    c("qunif", list(min=0.1, max=0.9), "abx.r"),           #probability of clearing R to become r
    c("qunif", list(min=0.1, max=0.9), "abx.s"),           #probability of clearing S to become s
    c("qunif", list(min=0, max=0.2), "repop.r"),           #probability of regrowth of s (sr—> sR)
    c("qunif", list(min=0, max=0.2), "repop.s1"),          #probability of regrowth of S  (s—>S)
    c("qunif", list(min=0, max=0.2), "repop.s2"),          #probability of regrowth of S  (sr—>Sr)
    c("qunif", list(min=0, max=0.2), "depop.r"),           #probability of sR-->sr without antibiotics
    c("qunif", list(min=3, max=7), "short_dur.s"),         #mean short duration of narrow spectrum antibiotics (normal distribution) 
    c("qunif", list(min=8, max=20), "long_dur.s"),         #mean long duration of narrow spectrum antibiotics (normal distribution) 
    c("qunif", list(min=3, max=7), "short_dur.r"),         #mean short duration of broad spectrum antibiotics (normal distribution) 
    c("qunif", list(min=8, max=20), "long_dur.r"),         #mean long duration of broad spectrum antibiotics (normal distribution) 
    c("qunif", list(min=1, max=2), "sdDur")                #standard deviation of the duration of antibiotics
    )

# arrange parameters in a way LHS will be happy with
q <- unlist(lapply(parameters, function(l) l[[1]]))
q.arg <- lapply(parameters, function(l) l[2:3])
factors <- unlist(lapply(parameters, function(l) l[[4]]))

# Test
source(paste0("model_binary.R"))
# if they don't follow the exact listing of function variables, they seem to feed the wrong range to the wrong variable...
# MAKE SURE the variable listing and ORDER MATCHES the variable listing input into diff_prevalence
if(!(sum(factors == parameters_binary) ==  length(parameters_binary))){
    stop("Test Error: Listing of parameters in cobweb does not match parameters accepted by diff_prevalence function.")
}

old <- Sys.time() # get start time
N=1000
LHS.binary <- LHS(modelRun.binary, factors, N=N, q, q.arg, nboot=1000, cl=cl)
results.binary <- get.results(LHS.binary)
new <- Sys.time() - old # calculate difference
print(new) # print in nice format

old <- Sys.time() # get start time
LHS.binary2 <- LHS(modelRun.binary, factors, N=N-100, q, q.arg, nboot=1000, cl=cl)
results.binary2 <- get.results(LHS.binary2)
# print elapsed time
new <- Sys.time() - old # calculate difference
print(new) # print in nice format

# Save run to disk
image_name <- paste0("./runs/LHS_", model, "_", N, format(Sys.time(), "%d%b%Y_%H%M%Z"))
save(LHS.binary, file=paste0(image_name, ".Rdata"))

image_name <- paste0("LHS2_", model, "_", N-100, format(Sys.time(), "%d%b%Y_%H%M%Z"))
save(LHS.binary2, file=paste0("./runs/", image_name, ".Rdata"))

stopCluster(cl)

