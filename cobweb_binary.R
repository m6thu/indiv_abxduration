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
    "mean.max.los",       #mean of length of stay (normal distribution)
    "p.s",                #probability of being prescribed narrow spectrum antibiotic
    "p.r.day1",           #probability of being prescribed broad spectrum antibiotic on day 1 of admission 
    "p.r.dayafter",       #probability of being prescribed broad spectrum antibiotic after admission (daily probability)
    "prob_StartBact_R",   #probability of initial carriage of resistant organisms
    "prop_S_nonR",        #proportion of S in (S+s): prob_start_S <- prop_S_nonR*(1-prob_StartBact_R)
    "prop_Sr_inR",        #proportion of Sr in (r+R): prob_start_Sr <- prop_Sr_inR*prob_StartBact_R
    "prop_sr_inR",        #proportion of sr in (r+r): prob_start_sr <- prop_sr_inR*prob_StartBact_R
    "pi_r1",              #probability of being transmitted r to S (S—> Sr)
    "bif",                #bacterial interference factor (pi_r2 = pi_r1 * bif )
    "mu1",                #probability of being decolonised to S (Sr—> S) 
    "mu2",                #probability of being decolonised to S (sr—> s) 
    "abx.r",              #probability of clearing R to become r
    "abx.s",              #probability of clearing S to become s
    "repop.r1",           #probability of transmission of r to S (s—> Sr) 
    "repop.r2",           #probability of regrowth of s (sR—> sr)
    "repop.s1",           #probability of regrowth of S  (s—>S)
    "repop.s2",           #probability of regrowth of S  (sr—>Sr)
    "depop.r",            #probability of sR-->sr without antibiotics
    "short_dur.s",        #mean short duration of narrow spectrum antibiotics (normal distribution) 
    "long_dur.s",         #mean long duration of narrow spectrum antibiotics (normal distribution) 
    "short_dur.r",        #mean short duration of broad spectrum antibiotics (normal distribution) 
    "long_dur.r",         #mean long duration of broad spectrum antibiotics (normal distribution) 
    "sdDur"               #standard deviation of the duration of antibiotics
)  
q <- c(                   #distributions of parameters 
    "qunif",              #number of beds in the ward
    "qunif",              #mean of length of stay (normal distribution)
    "qunif",              #probability of being prescribed narrow spectrum antibiotic
    "qunif",              #probability of being prescribed broad spectrum antibiotic on day 1 of admission 
    "qunif",              #probability of being prescribed broad spectrum antibiotic after admission (daily probability)
    "qunif",              #probability of initial carriage of resistant organisms
    "qunif",              #proportion of S in (S+s): prob_start_S <- prop_S_nonR*(1-prob_StartBact_R)
    "qunif",              #proportion of Sr in (r+R): prob_start_Sr <- prop_Sr_inR*prob_StartBact_R
    "qunif",              #proportion of sr in (r+r): prob_start_sr <- prop_sr_inR*prob_StartBact_R
    "qunif",              #probability of being transmitted r to S (S—> Sr)
    "qunif",              #bacterial interference factor (pi_r2 = pi_r1 * bif )
    "qunif",              #probability of being decolonised to S (Sr—> S) 
    "qunif",              #probability of being decolonised to S (sr—> s) 
    "qunif",              #probability of clearing R to become r
    "qunif",              #probability of clearing S to become s
    "qunif",              #probability of transmission of r to S (s—> Sr) 
    "qunif",              #probability of regrowth of s (sR—> sr)
    "qunif",              #probability of regrowth of S  (s—>S)
    "qunif",              #probability of regrowth of S  (sr—>Sr)
    "qunif",              #probability of probability of sR-->sr without antibiotics
    "qunif",              #mean short duration of narrow spectrum antibiotics (normal distribution) 
    "qunif",              #mean long duration of narrow spectrum antibiotics (normal distribution) 
    "qunif",              #mean short duration of broad spectrum antibiotics (normal distribution) 
    "qunif",              #mean long duration of broad spectrum antibiotics (normal distribution) 
    "qunif"               #standard deviation of the duration of antibiotics
) 

q.arg <- list(                  #set limits of parameters 
    list(min=5, max=50),        #n.bed; number of beds in the ward
    list(min=3, max=10),        #mean.max.los; mean of length of stay (normal distribution)
    list(min=0.01, max=0.5),    #probability of being prescribed narrow spectrum antibiotic
    list(min=0.01, max=0.5),    #probability of being prescribed broad spectrum antibiotic on day 1 of admission 
    list(min=0.00001, max=0.5), #probability of being prescribed broad spectrum antibiotic after admission (daily probability)
    list(min=0.01, max=0.99),   #probability of initial carriage of resistant organisms
    list(min=0.01, max=0.99),   #proportion of S in (S+s): prob_start_S <- prop_S_nonR*(1-prob_StartBact_R)
    list(min=0.01, max=0.99),   #proportion of Sr in (r+R): prob_start_Sr <- prop_Sr_inR*prob_StartBact_R
    list(min=0.01, max=0.99),   #proportion of sr in (r+r): prob_start_sr <- prop_sr_inR*prob_StartBact_R
    list(min=0.0001, max=0.001), #pi_r1; probability of being transmitted r to S (S—> Sr)
    list(min=0, max=0.2),       #bif; bacterial interference factor (pi_r2 = pi_r1 * bif )
    list(min=0, max=0.03),      #probability of being decolonised to S (Sr—> S) 
    list(min=0, max=0.03),      #probability of being decolonised to S (sr—> s) 
    list(min=0.01, max=0.9),    #probability of clearing R to become r
    list(min=0.01, max=0.9),    #probability of clearing S to become s
    list(min=0.01, max=0.9),    #probability of transmission of r to S (s—> Sr) 
    list(min=0.0001, max=20),   #probability of sR—> sr
    list(min=0, max=0.03),      #repop.s1; probability of regrowth of S (s—>S)
    list(min=0, max=0.03),      #probability of regrowth of S (sr—>Sr)
    list(min=0, max=0.03),      #probability of sR-->sr without antibiotics
    list(min=3, max=7),         #mean short duration of narrow spectrum antibiotics (normal distribution) 
    list(min=3, max=7),         #mean long duration of narrow spectrum antibiotics (normal distribution) 
    list(min=10, max=21),       #mean short duration of broad spectrum antibiotics (normal distribution) 
    list(min=10, max=21),       #mean long duration of broad spectrum antibiotics (normal distribution) 
    list(min=1, max=5 )         #standard deviation of the duration of antibiotics 
)

modelRun.binary <- function (data.df) { #data.df is a dataframe of the parameter values in columns 
    return(mapply(diff_prevalence, 
                  data.df[,1], data.df[,2], data.df[,3], data.df[,4], data.df[,5], 
                  data.df[,6], data.df[,7], data.df[,8], data.df[,9], 
                  data.df[,10], data.df[,11], data.df[,12], data.df[,13], data.df[,14], data.df[,15], 
                  data.df[,16], data.df[,17], data.df[,18],data.df[,19],data.df[,20], 
                  data.df[,21], data.df[,22], data.df[,23], data.df[,24], data.df[,25]
                  
    ))
}

# Use the LHD function to generate a hypercube
old <- Sys.time() # get start time
LHS.binary<- LHS(modelRun.binary, factors, 40, q, q.arg, nboot=10)
# print elapsed time
new <- Sys.time() - old # calculate difference
print(new) # print in nice format
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
require(plotrix) #load MASS package
yellow<-alpha("#E69F00", alpha=0.3)
green<-alpha("#009E73", alpha=0.6)
colors<- c(yellow, green) #choose 2 colors - 1 for parameters that produced top 5% of outcomes and one for the rest
outcome.df$top5<- as.factor(outcome.df$top5)
parcoordlabel<-function (x, col = 1, lty = 1,  lblcol="black",...) 
{
    df <- as.data.frame(x)
    pr <- lapply(df, pretty)
    rx <- lapply(pr, range, na.rm = TRUE)
    x <- mapply(function(x,r) {
        (x-r[1])/(r[2]-r[1])
    },
    df, rx)
    matplot(1L:ncol(x), t(x), type = "l", col = col, lty = lty, 
            xlab = "", ylab = "", axes = FALSE,...)
    axis(1, at = 1L:ncol(x), labels = c(colnames(x)), las = 2)
    for (i in 1L:ncol(x)) {
        lines(c(i, i), c(0, 1), col = "grey")
        text(c(i, i), seq(0,1,length.out=length(pr[[i]])), labels = pr[[i]], 
             xpd = NA, col=lblcol, cex=0.5)
    }
    invisible()
}
parcoordlabel(outcome.df[,c(1:length(factors))], col = colors[outcome.df$top5])

#5. Check agreement between runs to decide if our sample size for adequate 
# Symmetric Blest Measure of Agreement (SBMA) between the PRCC coeffients of two runs with different sample sizes.
check.LHS.binary <- LHS(modelRun.binary, factors.binary, 100, q.binary, q.arg.binary)
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