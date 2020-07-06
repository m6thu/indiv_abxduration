###################################################################################
###Effect of antibiotic duration on resistance carriage in hospitalised patients###
###### explore magnitude of difference in carrier prevalence in simple model ######
###################################################################################

# SAMPLE PARAMETER SPACE 
# load libraries 
require(pse) #load pse package for Latin Hypercube
require(sensitivity) #load sensitivity package for sensitivity analysis 
require(parallel) # load parallel processing package to use multiple cores on computer (or cluster)

rm(list=ls()) # Clean working environment

cl <- makeCluster(detectCores(), outfile=paste0('/Users/moyin/Desktop/parallel_error_simple.txt'))

model <- 'simple'

# source functions on all cores
clusterCall(cl, function() {source('model_simple.R')})

modelRun.simple <- function (data.df) { #data.df is a dataframe of the parameter values in columns 
    return(mapply(prevalence, 
                  data.df[,1], data.df[,2], data.df[,3], 
                  data.df[,4], data.df[,5], data.df[,6], 
                  data.df[,7], data.df[,8], data.df[,9], 
                  data.df[,10], data.df[,11], data.df[,12], 
                  data.df[,13],  data.df[,14]
    ))
}

################################## Define parameters and run LHS ####################################
#list parameters, the probability density functions from which the parameter values will be calculated, and what are the arguments to these functions
#list parameters together with name, so they are "linked" or not easily confused
parameters <- list(
  c("qunif", list(min=15, max=15.01), "n.bed"),        #"n.bed", number of beds in the ward
  c("qunif", list(min=3, max=7), "max.los"),           #"max.los", mean of length of stay
  c("qunif", list(min=0.6, max=0.8), "prop_R"),        #"prob_StartBact_R",probability of initial carriage of resistant organisms
  c("qunif", list(min=0, max=1), "prop_S"),            #"prop_S", proportion of S in the population of S and ss
  c("qunif", list(min=0, max=1), "bif"),               #"bif", bacterial interference factor
  c("qunif", list(min=0.3, max=0.301), "pi_ssr"),      # "pi_ssr" probability of being transmitted r to ss (ss—> ssr)
  c("qunif", list(min=0.02, max=0.12), "repop.s"),    # "repop.S" probability of ss repopulated to S (Palleja, Nature Biology, 2018 on gut recovery ~9 months)
  c("qunif", list(min=0.002, max=0.02), "mu"),         # "mu", probability of decolonisation (Haggai Bar-Yoseph, JAC, 2016, decreasing colonization rates from 76.7% (95% CI=69.3%–82.8%) at 1 month to 35.2% (95% CI=28.2%–42.9%) at 12 months of follow-up)
  c("qunif", list(min=0.1, max=0.5), "abx.s"),         # "abx.s", probability of S becoming ss after being on narrow spectrum antibiotics
  c("qunif", list(min=0, max=0.0000000001), "abx.r"),  # "abx.r", probability of R becoming ss after being on broad spectrum antibiotics
  c("qunif", list(min=0.8, max=0.801), "p.infect"),    # "p.infect", probability of being prescribed antibiotics
  c("qunif", list(min=30, max=100), "cum.r.1"),        # admission day when cummulative prabability of HAI requiring abx.r is 1
  c("qunif", list(min=0.8, max=0.801), "p.r.day1"),    # probability of being prescribed broad spectrum antibiotic on admission 
  c("qunif", list(min=3, max=21), "meanDur")           # "meanDur", mean duration of antibiotics (normal distribution)
)  

# arrange parameters in a way LHS will be happy with
q <- unlist(lapply(parameters, function(l) l[[1]]))
q.arg <- lapply(parameters, function(l) l[2:3])
factors <- unlist(lapply(parameters, function(l) l[[4]]))

# Test
# if they don't follow the exact listing of function variables, they seem to feed the wrong range to the wrong variable...
# MAKE SURE the variable listing and ORDER MATCHES the variable listing input into diff_prevalence
source(paste0("model_simple.R"))
if(!(sum(factors == parameters_prevalence_simple) ==  length(parameters_prevalence_simple))){
    stop("Test Error: Listing of parameters in cobweb does not match parameters accepted by prevalence function.")
}

# Use the LHD function to generate a hypercube 
## run 1 
abxr='zero_high'
old <- Sys.time() # get start time
N = 50
LHS.simple <- LHS(modelRun.simple, factors, N=N, q, q.arg, nboot=100, cl=cl) #N is the size of the hypercube
# print elapsed time
new <- Sys.time() - old # calculate difference
print(new) # print in nice format
image_name <- paste0("LHS_", model, "_", N, "_",abxr,format(Sys.time(), "%d%b%Y_%H%M%Z"))
save(LHS.simple, file=paste0("./runs/", image_name, ".Rdata"))

## run 2 
N=N+100
old <- Sys.time() # get start time
LHS.simple2 <- LHS(modelRun.simple, factors, N=N, q, q.arg, nboot=100, cl=cl)
# print elapsed time
new <- Sys.time() - old # calculate difference
print(new) # print in nice format
image_name <- paste0("LHS_", model, "_", N, "_",abxr,format(Sys.time(), "%d%b%Y_%H%M%Z"))
save(LHS.simple2, file=paste0("./runs/", image_name, ".Rdata"))

stopCluster(cl)
