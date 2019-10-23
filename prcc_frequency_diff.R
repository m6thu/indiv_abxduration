################################### Dependencies and functions ################################################

# SAMPLE PARAMETER SPACE 
# load libraries 
require(pse) #load pse package for Latin Hypercube
require(sensitivity) #load sensitivity package for sensitivity analysis 
require(parallel) # load parallel processing package to use multiple cores on computer (or cluster)

setwd('/Users/moyin/Documents/nBox/git_projects/indiv_abxduration/')

cl <- makeCluster(detectCores())

model <- 'frequency'
# source functions on all cores
clusterCall(cl, function() {source('model_frequency.R')})

modelRun.freq <- function (data.df) { #data.df is a dataframe of the parameter values in columns 
    return(mapply(diff_prevalence, 
                  data.df[,1], data.df[,2], data.df[,3], 
                  data.df[,4], data.df[,5], data.df[,6], 
                  data.df[,7], data.df[,8], data.df[,9], 
                  data.df[,10], data.df[,11], data.df[,12], 
                  data.df[,13], data.df[,14], data.df[,15], 
                  data.df[,16], data.df[,17]
    ))
}

################################## Define parameters and run LHS ####################################
#list parameters, the probability density functions from which the parameter values will be calculated, and what are the arguments to these functions
# Also note that pse gets fussy if min=max and will complain
#Error in if (max(abs(cor(vars)[l, 1:(l - 1)] - COR[l, 1:(l - 1)])) < eps) { : 
#        missing value where TRUE/FALSE needed
parameters <- list(
    c("qunif", list(min=5, max=50), "n.bed"),              # n.bed; number of beds in the ward
    c("qunif", list(min=3, max=20), "max.los"),       # max.los; mean of length of stay (normal distribution)
    c("qunif", list(min=0.1, max=1), "p.infect"),          # probability of being prescribed narrow spectrum antibiotic
    c("qunif", list(min=10, max=1000), "cum.r.1"),        # admission day when cummulative prabability of HAI requiring abx.r is 1
    c("qunif", list(min=0.1, max=1), "p.r.day1"),          # probability of being prescribed broad spectrum antibiotic on day 1 of admission 
    c("qunif", list(min=2, max=27), "K"),                  # gut holding capacity, on log scale, largest R number possible is exp(300) - typical colonic bacteria 10^14 number/mL content https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4991899/
    c("qunif", list(min=0.1, max=1), "total_prop"),      # mean of total starting amount of enterobacteriaceae on log scale
    c("qunif", list(min=0, max=1), "prop_R"),            # probability of a patient coming into the ward carrying R
    c("qunif", list(min=0,max=0.002), "pi_ssr"),             # pi_ssr = daily probability of transmitting resistant E coli
    c("qunif", list(min=0,max=1), "r_mean"),             # r_mean = mean R proportion for those carrying R
    c("qunif", list(min=0.01,max=0.05), "r_growth"),         # r_growth = growth constant for logistic growth
    c("qunif", list(min=0,max=0.9), "r_thres"),             # r_thres = threshold amount of bacteria before R can be transmitted
    c("qunif", list(min=0.005,max=0.015), "s_growth"),         # s_growth = amount transmitted on log scale
    c("qunif", list(min=10,max=15), "abx.s"),               # abxr_killr = amount of r killed by broad spectrum abx r
    c("qunif", list(min=10,max=15), "abx.r"),               # abxr_kills = amount of s killed by broad spectrum abx r
    c("qunif", list(min=3, max=7), "short_dur"),           # mean short duration of narrow spectrum antibiotics (normal distribution) 
    c("qunif", list(min=14, max=21), "long_dur")           # mean long duration of narrow spectrum antibiotics (normal distribution)
    )

# arrange parameters in a way LHS will be happy with
q <- unlist(lapply(parameters, function(l) l[[1]]))
q.arg <- lapply(parameters, function(l) l[2:3])
factors <- unlist(lapply(parameters, function(l) l[[4]]))

# Test
# if they don't follow the exact listing of function variables, they seem to feed the wrong range to the wrong variable...
# MAKE SURE the variable listing and ORDER MATCHES the variable listing input into diff_prevalence
source(paste0("model_frequency.R"))
if(!(sum(factors == parameters_diff_prevalence_freq) ==  length(parameters_diff_prevalence_freq))){
    stop("Test Error: Listing of parameters in cobweb does not match parameters accepted by diff_prevalence function.")
}

# Use the LHS function to generate a hypercube 
##run 1
abxr='notzero'
N=1000
old <- Sys.time() # get start time
LHS.freq<- LHS(modelRun.freq, factors, N=N, q, q.arg, nboot=100,cl=cl)
# print elapsed time
new <- Sys.time() - old # calculate difference
print(new) # print in nice format
# Save run to disk
image_name <- paste0("LHSdiff_", model, "_", N, "_",abxr,"_",format(Sys.time(), "%d%b%Y_%H%M%Z"))
save(LHS.freq,file=paste0("./runs/", image_name, ".Rdata"))

##run 2
N=N+100
old <- Sys.time() # get start time
LHS.freq2 <- LHS(modelRun.freq, factors, N=N, q, q.arg, nboot=100, cl=cl)
new <- Sys.time() - old # calculate difference
print(new) # print in nice format
# Save run to disk
image_name <- paste0("LHSdiff_", model, "_", N,"_",abxr,"_",format(Sys.time(), "%d%b%Y_%H%M%Z"))
save(LHS.freq2,file=paste0("./runs/", image_name, ".Rdata"))

stopCluster(cl)
