#######Modelling Day Project######
#######Parameter exploration######
################################### Dependencies and functions ################################################
setwd('/Users/moyin/Documents/nBox/git_projects/indiv_abxduration/')

# SAMPLE PARAMETER SPACE 
# load libraries 
require(pse) #load pse package for Latin Hypercube
require(sensitivity) #load sensitivity package for sensitivity analysis
require(parallel) # load parallel processing package to use multiple cores on computer (or cluster)

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
                  data.df[,16], data.df[,17], data.df[,18]
    ))
}

################################## Define parameters and run LHS ####################################
#list parameters, the probability density functions from which the parameter values will be calculated, and what are the arguments to these functions
#list parameters together with name, so they are "linked" or not easily confused
parameters <- list(
    c("qunif", list(min=5, max=50), "n.bed"),              #n.bed; number of beds in the ward
    c("qunif", list(min=3, max=20), "max.los"),       #max.los; mean of length of stay (exponential distribution)
    c("qunif", list(min=0, max=1), "prop_R"),    #probability of initial carriage of resistant organisms
    c("qunif", list(min=0, max=1), "prop_S_nonR"),         #proportion of S in (S+s): prob_start_S <- prop_S_nonR*(1-prob_R)
    c("qunif", list(min=0, max=1), "prop_Sr_inR"),         #proportion of Sr in (r+R): prob_start_Sr <- prop_Sr_inR*prob_R
    c("qunif", list(min=0, max=1), "prop_sr_inR"),         #proportion of sr in (r+r): prob_start_sr <- prop_sr_inR*prob_R
    c("qunif", list(min=0, max=1), "bif"),                 #bacterial interference factor (pi_ssr = pi_r1 * bif )
    c("qunif", list(min=0, max=0.002), "pi_ssr"),            #probability of being transmitted r to ss (ss—> ssr)
    c("qunif", list(min=0.002, max=0.02), "repop.s"),     #probability of regrowth of S  (s—>S)
    c("qunif", list(min=0.01, max=0.05), "repop.r"),     #probability of regrowth of s (sr—> sR)
    c("qunif", list(min=0.002, max=0.02), "mu"),          #probability of being decolonised to S (Sr—> S) 
    c("qunif", list(min=0.1, max=0.5), "abx.s"),           #probability of clearing S to become s
    c("qunif", list(min=0, max=0.00000001), "abx.r"),           #probability of clearing R to become r
    c("qunif", list(min=0.1, max=1), "p.infect"),          #probability of being prescribed narrow spectrum antibiotic
    c("qunif", list(min=10, max=1000), "cum.r.1"),        #admission day when cummulative prabability of HAI requiring abx.r is 1
    c("qunif", list(min=0.1, max=1), "p.r.day1"),          #probability of being prescribed broad spectrum antibiotic on day 1 of admission 
    c("qunif", list(min=3, max=7), "short_dur"),           #mean short duration of antibiotics (normal distribution) 
    c("qunif", list(min=14, max=21), "long_dur")           #mean long duration of antibiotics (normal distribution) 
    )

# arrange parameters in a way LHS will be happy with
q <- unlist(lapply(parameters, function(l) l[[1]]))
q.arg <- lapply(parameters, function(l) l[2:3])
factors <- unlist(lapply(parameters, function(l) l[[4]]))

# Test
source(paste0("model_binary.R"))
# if they don't follow the exact listing of function variables, they seem to feed the wrong range to the wrong variable...
# MAKE SURE the variable listing and ORDER MATCHES the variable listing input into diff_prevalence
if(!(sum(factors == parameters_diff_prevalence_binary) ==  length(parameters_diff_prevalence_binary))){
    stop("Test Error: Listing of parameters in cobweb does not match parameters accepted by diff_prevalence function.")
}


#Run model and save runs 
##run 1
abxr='zero'
old <- Sys.time() # get start time
N=1200
LHS.binary <- LHS(modelRun.binary, factors, N=N, q, q.arg, nboot=100, cl=cl)
new <- Sys.time() - old # calculate difference
print(new) # print in nice format
# Save run to disk
image_name <- paste0("LHSdiff_", model, "_", N, "_",abxr, "_",format(Sys.time(), "%d%b%Y_%H%M%Z"))
save(LHS.binary, file=paste0("./runs/", image_name, ".Rdata"))

N=N+100
old <- Sys.time() # get start time
LHS.binary2 <- LHS(modelRun.binary, factors, N=N, q, q.arg, nboot=100, cl=cl)
# print elapsed time
new <- Sys.time() - old # calculate difference
print(new) # print in nice format
image_name <- paste0("LHSdiff_", model, "_", N, "_",abxr, "_",format(Sys.time(), "%d%b%Y_%H%M%Z"))
save(LHS.binary2, file=paste0("./runs/", image_name, ".Rdata"))

stopCluster(cl)

