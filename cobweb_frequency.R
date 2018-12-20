#######Modelling Day Project######
#######Parameter exploration######
################################### Dependencies and functions ################################################

# SAMPLE PARAMETER SPACE 
# load libraries 
require(pse) #load pse package for Latin Hypercube
require(sensitivity) #load sensitivity package for sensitivity analysis 
require(parallel) # load parallel processing package to use multiple cores on computer (or cluster)

cl <- makeCluster(detectCores())

model <- 'frequency'
#source(paste0("model_frequency.R"))
# source functions on all cores
clusterCall(cl, function() {source('~/Desktop/indivi_duration/model_frequency.R')})

modelRun.freq <- function (data.df) { #data.df is a dataframe of the parameter values in columns 
    return(mapply(diff_prevalence, 
                  data.df[,1], data.df[,2], data.df[,3], 
                  data.df[,4], data.df[,5], data.df[,6], 
                  data.df[,7], data.df[,8], data.df[,9], 
                  data.df[,10], data.df[,11], data.df[,12], 
                  data.df[,13], data.df[,14], data.df[,15], 
                  data.df[,16], data.df[,17], data.df[,18], 
                  data.df[,19], data.df[,20], data.df[,21], data.df[,22]
    ))
}

################################## Define parameters and run LHS ####################################
#list parameters, the probability density functions from which the parameter values will be calculated, and what are the arguments to these functions
# Also note that pse gets fussy if min=max and will complain
#Error in if (max(abs(cor(vars)[l, 1:(l - 1)] - COR[l, 1:(l - 1)])) < eps) { : 
#        missing value where TRUE/FALSE needed
parameters <- list(
    c("qunif", list(min=5, max=50), "n.bed"), #n.bed; number of beds in the ward
    c("qunif", list(min=3, max=10), "mean.max.los"), #mean.max.los; mean of length of stay (normal distribution)
    c("qunif", list(min=0.01, max=0.5), "p.s"), #probability of being prescribed narrow spectrum antibiotic
    c("qunif", list(min=0.01, max=0.5), "p.r.day1"),           #probability of being prescribed broad spectrum antibiotic on day 1 of admission 
    c("qunif", list(min=0.00001, max=0.5), "p.r.dayafter"),       #probability of being prescribed broad spectrum antibiotic after admission (daily probability)
    c("qunif", list(min=100, max=10000), "K"), # gut holding capacity
    c("qunif", list(min=3, max=5), "t_mean"), # mean of total starting amount of gut bacteria on log scale
    c("qunif", list(min=0.5, max=2), "t_sd"),  # sd of total starting amount of gut bacteria on log scale
    c("qunif", list(min=1,max=2), "r_mean"), # mean of starting amount of resistant gut bacteria on log scale
    c("qunif", list(min=1,max=2), "r_sd"), #sd of starting amount of resistant gut bacteria on log scale
    c("qunif", list(min=0.01,max=0.10), "pi_r"), # pi_r = daily probability of transmitting resistant E coli
    c("qunif", list(min=5,max=20), "r_thres"), # r_thres = R threshold level for tranmissibility
    c("qunif", list(min=1,max=5), "r_growth"), # r_growth = growth constant for logistic growth
    c("qunif", list(min=5,max=20), "r_trans"), # r_trans = amount transmitted
    c("qunif", list(min=5,max=20), "abxr_killr"), # abxr_killr = amount of r killed by broad spectrum abx r
    c("qunif", list(min=5,max=20), "abxr_kills"), # abxr_kills = amount of s killed by broad spectrum abx r
    c("qunif", list(min=5,max=20), "abxs_kills"), # abxs_kills = amount of s killed by narrow spectrum abx s
    c("qunif", list(min=3, max=7), "short_dur.s"),        #mean short duration of narrow spectrum antibiotics (normal distribution) 
    c("qunif", list(min=3, max=7), "long_dur.s"),         #mean long duration of narrow spectrum antibiotics (normal distribution) 
    c("qunif", list(min=10, max=21), "short_dur.r"),        #mean short duration of broad spectrum antibiotics (normal distribution) 
    c("qunif", list(min=10, max=21), "long_dur.r"),         #mean long duration of broad spectrum antibiotics (normal distribution) 
    c("qunif", list(min=1, max=5), "sdDur")               #standard deviation of the duration of antibiotics
    )

# arrange parameters in a way LHS will be happy with
q <- unlist(lapply(parameters, function(l) l[[1]]))
q.arg <- lapply(parameters, function(l) l[2:3])
factors <- unlist(lapply(parameters, function(l) l[[4]]))

# Test
# if they don't follow the exact listing of function variables, they seem to feed the wrong range to the wrong variable...
# MAKE SURE the variable listing and ORDER MATCHES the variable listing input into diff_prevalence
if(!(sum(factors == parameters_frequency) ==  length(parameters_frequency))){
    stop("Test Error: Listing of parameters in cobweb does not match parameters accepted by diff_prevalence function.")
}

# Use the LHD function to generate a hypercube 
old <- Sys.time() # get start time
LHS.freq<- LHS(modelRun.freq, factors, 1000, q, q.arg, nboot=20)
# print elapsed time
new <- Sys.time() - old # calculate difference
print(new) # print in nice format


old <- Sys.time() # get start time
check.LHS.freq <- LHS(modelRun.freq, factors, 500, q, q.arg, nboot=10, cl=cl)
new <- Sys.time() - old # calculate difference
print(new) # print in nice format


results.freq <- get.results(LHS.freq)
# Save run to disk
image_name <- paste0("./runs/LHS_", model, "_", format(Sys.time(), "%d%b%Y_%H%M%Z"))
save.image(paste0(image_name, ".Rdata"))

stopCluster(cl)

##################################### Display results ########################################
#Plot findings: 
#1. empirical cumulative distribution function used to illustrate the distribution of the model results
plotecdf(LHS.freq) #outcome has a high probability in the steepest parts of the graph 

#2. scatterplot of the result as a function of each parameter: distribution of values returned by the model 
# in the parameter space sampled by the hypercube and how sensible are these model responses to the variation of each parameter.
plotscatter(LHS.freq)

#3. partial correlation coefficient measures how strong are the inear associations between the result 
# and each input parameter, after removing the linear effect of the other parameters. (CI generated by boostrapping)
plotprcc(LHS.freq)

pic(LHS.freq, nboot=40) #pic (partial inclination coefficient) is the sensitivity" of the model response in respect to each parameter
#represent the beta terms in y = alpha + beta*x regressions, after removing the linear effect of the other parameters

# 4. Cobweb
outcome.df<-as.data.frame(cbind(LHS.freq$data,results.freq)) #dummy matrix with parameter values in columns and outcome in last column
names(outcome.df)<- c(factors,'outcome') #name the columns of the dummy matrix 
for (i in 1:nrow(outcome.df)) {       #label the rows of parameter values that produced top 5% of the outcomes
    if (outcome.df$outcome[i]<quantile(outcome.df$outcome,probs = 0.95)) { 
        outcome.df$top5[i] <-0 } else {
            outcome.df$top5[i] <-1
        }
}
require(MASS) #load MASS package
colors<- c("#E69F00", "#009E73") #choose 2 colors - 1 for parameters that produced top 5% of outcomes and one for the rest
outcome.df$top5<- as.factor(outcome.df$top5)
parcoord(outcome.df[,c(1:length(factors))] , col= colors[outcome.df$top5],var.label=T)

#5. Check agreement between runs to decide if our sample size for adequate 
# Symmetric Blest Measure of Agreement (SBMA) between the PRCC coeffients of two runs with different sample sizes.
check.LHS.freq <- LHS(modelRun.freq, factors.freq, 250, q.freq, q.arg.freq)
(mySbma <- sbma(LHS.freq, check.LHS.freq))
# value of -1 indicates complete disagreement between the runs 
# value of 1 indicated complete agreement  (>0.7 acceptable)
#caveat: if none of the model parameters is monotonically correlated with the output, the agreement between runs may stay as low as 0.2 even for very large hypercubes.

#################################################################################
###DUMMY CODE FOR PARAMETER EXPLORATION###
fun.trial<- function(x1, x2, x3) {x1+x2^2+ x3^3} #dummy model 
factors.trial<-c('x1', 'x2', 'x3') #dummy variables 
q.trial<- c("qunif","qunif","qunif") #distribution of the variables 
q.arg.trial <- list(        #set limits to the variables   
    list(min=0.1, max=0.9),     
    list(min=0.1, max=0.9),       
    list(min=0.1, max=0.9)) 
modelRun.trial <- function (data.df) { #data.df is a dataframe of the parameter values in columns 
    return(mapply(fun.trial, data.df[,1], data.df[,2], data.df[,3]))
}

LHS.trial<- LHS(modelRun.trial, factors.trial, 50, q.trial, q.arg.trial, nboot=20) #50 parameter combinations to be generated, nboot= number of correlation coefficients 
results.trial<-get.results(LHS.trial)

#Plot findings: 
#1. empirical cumulative distribution function used to illustrate the distribution of the model results
plotecdf(LHS.trial) #outcome has a high probability in the steepest parts of the graph 

#2. scatterplot of the result as a function of each parameter: distribution of values returned by the model 
# in the parameter space sampled by the hypercube and how sensible are these model responses to the variation of each parameter.
plotscatter(LHS.trial)

#3. partial correlation coefficient measures how strong are the inear associations between the result 
# and each input parameter, after removing the linear effect of the other parameters. (CI generated by boostrapping)
plotprcc(LHS.trial)

pic(LHS.trial, nboot=40) #pic (partial inclination coefficient) is the sensitivity" of the model response in respect to each parameter
#represent the beta terms in y = alpha + beta*x regressions, after removing the linear effect of the other parameters

# 4. Cobweb
outcome.df<-as.data.frame(cbind(LHS.trial$data,results.trial)) #dummy matrix with parameter values in columns and outcome in last column
names(outcome.df)<- c('X1','X2','X3','Y') #name the columns of the dummy matrix 
for (i in 1:nrow(outcome.df)) {       #label the rows of parameter values that produced top 5% of the outcomes
    if (outcome.df$Y[i]<quantile(outcome.df$Y,probs = 0.95)) { 
        outcome.df$top5[i] <-0 } else {
            outcome.df$top5[i] <-1
        }
}
library(MASS) #load MASS package
colors<- c("#E69F00", "#009E73") #choose 2 colors - 1 for parameters that produced top 5% of outcomes and one for the rest
outcome.df$top5<- as.factor(outcome.df$top5)
parcoord(outcome.df[,c(1:3)] , col= colors[outcome.df$top5],var.label=T)

#5. Check agreement between runs to decide if our sample size for adequate 
# Symmetric Blest Measure of Agreement (SBMA) between the PRCC coeffients of two runs with different sample sizes.
#check.LHS.trial <- LHS(modelRun.trial, factors.trial, 250, q.trial, q.arg.trial)
(mySbma <- sbma(LHS.trial, check.LHS.trial))
# value of -1 indicates complete disagreement between the runs 
# value of 1 indicated complete agreement  (>0.7 acceptable)
#caveat: if none of the model parameters is monotonically correlated with the output, the agreement between runs may stay as low as 0.2 even for very large hypercubes.