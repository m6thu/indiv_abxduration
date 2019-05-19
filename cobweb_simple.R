################################### Dependencies and functions ################################################

# SAMPLE PARAMETER SPACE 
# load libraries 
require(pse) #load pse package for Latin Hypercube
require(sensitivity) #load sensitivity package for sensitivity analysis 
require(parallel) # load parallel processing package to use multiple cores on computer (or cluster)

setwd('/Users/moyin/Desktop/angelsfly/indiv_abxduration')

cl <- makeCluster(detectCores())

model <- 'simple'
# source functions on all cores
clusterCall(cl, function() {source('model_simple.R')})

modelRun.simple <- function (data.df) { #data.df is a dataframe of the parameter values in columns 
    return(mapply(diff_prevalence, 
                  data.df[,1], data.df[,2], data.df[,3], 
                  data.df[,4], data.df[,5], data.df[,6], 
                  data.df[,7], data.df[,8], data.df[,9], 
                  data.df[,10], data.df[,11], data.df[,12], 
                  data.df[,13]
    ))
}

################################## Define parameters and run LHS ####################################
#list parameters, the probability density functions from which the parameter values will be calculated, and what are the arguments to these functions
#list parameters together with name, so they are "linked" or not easily confused
parameters <- list(
    c("qunif", list(min=3, max=50), "n.bed"),  #"n.bed", number of beds in the ward
    c("qunif", list(min=3, max=30), "mean.max.los"), #"mean.max.los", mean of length of stay
    c("qunif", list(min=0, max=0.9), "prob_StartBact_R"),   #"prob_StartBact_R",probability of initial carriage of resistant organisms
    c("qunif", list(min=0, max=1), "prop_S_nonR"),        #"prop_S_nonR", proportion of S in the population of S and ss
    c("qunif", list(min=0, max=1), "bif"),                #"bif", bacterial interference factor
    c("qunif", list(min=0, max=0.05), "pi_ssr"),              # "pi_ssr" probability of being transmitted r to ss (ss—> ssr)
    c("qunif", list(min=0, max=0.2), "repop.s1"),          # "repop.s1" probability of ss repopulated to S (Palleja, Nature Biology, 2018 on gut recovery ~9 months)
    c("qunif", list(min=0, max=0.005), "mu_r"),                 # "mu_r", probability of decolonisation (Haggai Bar-Yoseph, JAC, 2016, decreasing colonization rates from 76.7% (95% CI=69.3%–82.8%) at 1 month to 35.2% (95% CI=28.2%–42.9%) at 12 months of follow-up)
    c("qunif", list(min=0.1, max=0.9), "abx.clear"),             # "abx.clear", probability of S becoming ss after being on antibiotics
    c("qunif", list(min=0.1, max=0.9), "p"),               # "p", probability of being prescribed antibiotics
    c("qunif", list(min=3, max=7), "short_dur"),           # "short_dur", mean short duration of antibiotics (normal distribution)
    c("qunif", list(min=8, max=20), "long_dur"),         # "long_dur", mean long duration of antibiotics (normal distribution)
    c("qunif", list(min=1, max=2), "sdDur")               # "sdDur", standard deviation of duration of antibiotics
)  

# arrange parameters in a way LHS will be happy with
q <- unlist(lapply(parameters, function(l) l[[1]]))
q.arg <- lapply(parameters, function(l) l[2:3])
factors <- unlist(lapply(parameters, function(l) l[[4]]))

# Test
# if they don't follow the exact listing of function variables, they seem to feed the wrong range to the wrong variable...
# MAKE SURE the variable listing and ORDER MATCHES the variable listing input into diff_prevalence
source(paste0("model_simple.R"))
if(!(sum(factors == parameters_simple) ==  length(parameters_simple))){
    stop("Test Error: Listing of parameters in cobweb does not match parameters accepted by diff_prevalence function.")
}

# Use the LHD function to generate a hypercube 
old <- Sys.time() # get start time
LHS.simple <- LHS(modelRun.simple, factors, N=20, q, q.arg, nboot=1000, cl=cl) #N is the size of the hypercube
# print elapsed time
new <- Sys.time() - old # calculate difference
print(new) # print in nice format

samples<-LHS.simple[['data']]
results <- get.results(LHS.simple)
aa_data_simple<-cbind(samples, results)

# Save run to disk
iterationstotry= c(1, 5, 50, 100, 300)
numberofrepeatsineachiteration=20
aa_data_simple<-list()
for (i in 1: (max(iterationstotry)*numberofrepeatsineachiteration)){
    LHS.simple <- LHS(modelRun.simple, factors, N=20, q, q.arg, nboot=1000, cl=cl)
    samples<-LHS.simple[['data']]
    results <- get.results(LHS.simple)
    aa_data_simple[[i]]<-cbind(samples, results)
} 

for (i in iterationstotry){
  data=total
  setwd('/Users/moyin/Desktop/angelsfly/indiv_abxduration/papers/simple_model')
  numbertoputinfolderi=i*numberofrepeatsineachiteration
  dir.create(as.character(i))
  for (k in 1:numberofrepeatsineachiteration){
    directory=paste0('/Users/moyin/Desktop/angelsfly/indiv_abxduration/papers/simple_model/',i)
    setwd(directory)
    dir.create(as.character(k))
    for (g in 1:i){
      directory=paste0('/Users/moyin/Desktop/angelsfly/indiv_abxduration/papers/simple_model/',i,'/',k)
      setwd(directory)
      dir.create(as.character(g))
      directory=paste0('/Users/moyin/Desktop/angelsfly/indiv_abxduration/papers/simple_model/',i,'/',k, '/',g)
      setwd(directory)
      write.csv(data[[1]],'aa_data_simple.csv',row.names=FALSE)
      data=data[-1]
    }
  }
}

image_name <- paste0("LHS_", model, "_", format(Sys.time(), "%d%b%Y_%H%M%Z"))
save.image(paste0("./runs/", image_name, ".Rdata"))

stopCluster(cl)
