#######Modelling Day Project######
#######Parameter exploration######
##################################

# SAMPLE PARAMETER SPACE 
# load libraries 
require(pse) #load pse package for Latin Hypercube
require(sensitivity) #load sensitivity package for sensitivity analysis 

# model can be "simple", "binary", or "frequency"
model <- "binary"

source(paste0("model_", model,".R"))

#list parameters, the probability density functions from which the parameter values will be calculated, and what are the arguments to these functions
factors <- c(             #list parameters in an array 
    "n.bed",              #number of beds in the ward
    #"n.days",             #number of days of observation 
    "mean.max.los",       #mean of length of stay (normal distribution)
    "p.s",                #probability of being prescribed narrow spectrum antibiotic
    "p.r",                #probability of being prescribed broad spectrum antibiotic
    "prob_StartBact_R",   #probability of initial carriage of resistant organisms
    "pi_r1",              #probability of being transmitted r to S (S—> Sr)
    "pi_r2",              #probability of being transmitted r to s (s—>sr)
    "mu1",                #probability of being decolonised to S (Sr—> S) 
    "mu2",                #probability of being decolonised to S (sr—> s) 
    "abx.r",              #probability of clearing R to become r
    "abx.s",              #probability of clearing S to become s
    "repop.r1",           #probability of transmission of r to S (s—> Sr) 
    "repop.r2",           #probability of regrowth of s (sR—> sr)
    "repop.s1",           #probability of regrowth of S  (s—>S)
    "repop.s2",           #probability of regrowth of S  (sr—>Sr)
    "repop.s3",           #probability of transmission of r to S (s—> Sr) 
    #"bif",               #bacterial interference factor (pi_r2 = pi_r1 * bif )
    #"bif1",              #bacterial interference factor 1 (repop.r3=repop.r2*bif1)
    "short_dur",          #mean short duration of antibiotics (normal distribution) 
    "long_dur"            #mean long duration of antibiotics (normal distribution) 
)  
q <- c(                   #distributions of parameters 
    "qunif",              # "n.bed", number of beds in the ward
    #"qunif",              # "n.days" number of days of observation 
    "qunif",              # "mean.max.los", mean of length of stay (normal distribution)
    "qunif",              # "p.s", probability of being prescribed narrow spectrum antibiotic
    "qunif",              # "p.r", probability of being prescribed broad spectrum antibiotic
    "qunif",              # "prob_StartBact_R",probability of initial carriage of resistant organisms
    "qunif",              # "pi_r1" probability of being transmitted r to S (S—> Sr)
    "qunif",              # "pi_r2",probability of being transmitted r to s (s—>sr)
    "qunif",              # "mu1" probability of being decolonised to S (Sr—> S) 
    "qunif",              # "mu2",probability of being decolonised to S (sr—> s) 
    "qunif",              # "abx.r", probability of clearing R to become r
    "qunif",              # "abx.s", probability of clearing S to become s
    "qunif",              # "repop.r1", probability of transmission of r to S (s—> Sr) 
    "qunif",              # "repop.r2", probability of regrowth of s (sR—> sr)
    "qunif",              # "repop.s1", probability of regrowth of S  (s—>S)
    "qunif",              # "repop.s2", probability of regrowth of S  (sr—>Sr)
    "qunif",              # "repop.s3",probability of transmission of r to S (s—> Sr) 
    "qunif",              #"short_dur", mean short duration of antibiotics (normal distribution) 
    "qunif"               #"long_dur", mean long duration of antibiotics (normal distribution)  
) 

q.arg <- list(            #set limits of parameters 
    list(min=10, max=20),                # "n.bed", number of beds in the ward
    #list(min=10, max=20),                # "n.days" number of days of observation 
    list(min=5, max=10),                 # "mean.max.los", mean of length of stay (normal distribution)
    list(min=0.1, max=0.9),              # "p.s", probability of being prescribed narrow spectrum antibiotic
    list(min=0.1, max=0.9),              # "p.r", probability of being prescribed broad spectrum antibiotic
    list(min=0.1, max=0.7),              # "prob_StartBact_R",probability of initial carriage of resistant organisms
    list(min=0.1, max=0.9),              # "pi_r1" probability of being transmitted r to S (S—> Sr)
    list(min=0.1, max=0.9),              # "pi_r2",probability of being transmitted r to s (s—>sr)
    list(min=0.1, max=0.9),              # "mu1" probability of being decolonised to S (Sr—> S) 
    list(min=0.1, max=0.9),              # "mu2",probability of being decolonised to S (sr—> s) 
    list(min=0.1, max=0.9),              # "abx.r", probability of clearing R to become r
    list(min=0.1, max=0.9),              # "abx.s", probability of clearing S to become s
    list(min=0.1, max=0.9),              # "repop.r1", probability of transmission of r to S (s—> Sr) 
    list(min=0.1, max=0.9),              # "repop.r2", probability of regrowth of s (sR—> sr)
    list(min=0.1, max=0.9),              # "repop.s1", probability of regrowth of S  (s—>S)
    list(min=0.1, max=0.9),              # "repop.s2", probability of regrowth of S  (sr—>Sr)
    list(min=0.1, max=0.9),              # "repop.s3",probability of transmission of r to S (s—> Sr) 
    list(min=4, max=10),                 #"short_dur", mean short duration of antibiotics (normal distribution) 
    list(min=14, max=20)                 #"long_dur", mean long duration of antibiotics (normal distribution)  
)

modelRun.binary <- function (data.df) { #data.df is a dataframe of the parameter values in columns 
    return(mapply(diff_prevalence, 
                  data.df[,1], data.df[,2], data.df[,3], 
                  data.df[,4], data.df[,5], data.df[,6], 
                  data.df[,7], data.df[,8], data.df[,9], 
                  data.df[,10], data.df[,11], data.df[,12], 
                  data.df[,13], data.df[,14], data.df[,15], 
                  data.df[,16], data.df[,17], data.df[,18]
    ))
}

# Use the LHD function to generate a hypercube 
LHS.binary<- LHS(modelRun.binary, factors, 20, q, q.arg, nboot=20)
results.binary<-get.results(LHS.binary)

# Save run to disk
image_name <- paste0("LHS_", model, "_", format(Sys.time(), "%d%b%Y_%H%M%Z"))
save.image(paste0(image_name, ".Rdata"))

#Plot findings: 
#1. empirical cumulative distribution function used to illustrate the distribution of the model results
plotecdf(LHS.binary) #outcome has a high probability in the steepest parts of the graph 

#2. scatterplot of the result as a function of each parameter: distribution of values returned by the model 
# in the parameter space sampled by the hypercube and how sensible are these model responses to the variation of each parameter.
plotscatter(LHS.binary)

#3. partial correlation coefficient measures how strong are the inear associations between the result 
# and each input parameter, after removing the linear effect of the other parameters. (CI generated by boostrapping)
plotprcc(LHS.binary)

pic(LHS.binary, nboot=40) #pic (partial inclination coefficient) is the sensitivity" of the model response in respect to each parameter
#represent the beta terms in y = alpha + beta*x regressions, after removing the linear effect of the other parameters

# 4. Cobweb
outcome.df<-as.data.frame(cbind(LHS.binary$data,results.binary)) #dummy matrix with parameter values in columns and outcome in last column
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
check.LHS.binary <- LHS(modelRun.binary, factors.binary, 250, q.binary, q.arg.binary)
(mySbma <- sbma(LHS.binary, check.LHS.binary))
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
check.LHS.trial <- LHS(modelRun.trial, factors.trial, 250, q.trial, q.arg.trial)
(mySbma <- sbma(LHS.trial, check.LHS.trial))
# value of -1 indicates complete disagreement between the runs 
# value of 1 indicated complete agreement  (>0.7 acceptable)
#caveat: if none of the model parameters is monotonically correlated with the output, the agreement between runs may stay as low as 0.2 even for very large hypercubes.