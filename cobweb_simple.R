#######Modelling Day Project######
#######Parameter exploration######
##################################
################################### Dependencies and functions ################################################

# SAMPLE PARAMETER SPACE 
# load libraries 
require(pse) #load pse package for Latin Hypercube
require(sensitivity) #load sensitivity package for sensitivity analysis 
require(parallel) # load parallel processing package to use multiple cores on computer (or cluster)

cl <- makeCluster(detectCores())

model <- 'simple'
# source functions on all cores
clusterCall(cl, function() {source('~/Desktop/indiv_abxduration/model_simple.R')})

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
    c("qunif", list(min=5, max=50), "n.bed"),  #"n.bed", number of beds in the ward
    c("qunif", list(min=3, max=10), "mean.max.los"), #"mean.max.los", mean of length of stay 
    c("qunif", list(min=0.01, max=0.99), "prob_StartBact_R"),   #"prob_StartBact_R",probability of initial carriage of resistant organisms
    c("qunif", list(min=0.01, max=0.99), "prop_S_nonR"),        #"prop_S_nonR", proportion of S in the population of S and ss 
    c("qunif", list(min=0, max=1), "bif"),                #"bif", bacterial interference factor 
    c("qunif", list(min=0.0001, max=0.05), "pi_ssr"),              # "pi_ssr" probability of being transmitted r to ss (ss—> ssr)
    c("qunif", list(min=0, max=0.03), "repop.s1"),          # "repop.s1" probability of ss repopulated to S (Palleja, Nature Biology, 2018 on gut recovery ~9 months)
    c("qunif", list(min=0, max=0.03), "mu_r"),                 # "mu_r", probability of decolonisation (Haggai Bar-Yoseph, JAC, 2016, decreasing colonization rates from 76.7% (95% CI = 69.3%–82.8%) at 1 month to 35.2% (95% CI = 28.2%–42.9%) at 12 months of follow-up)
    c("qunif", list(min=0.1, max=0.9), "abx.clear"),             # "abx.clear", probability of S becoming ss after being on antibiotics
    c("qunif", list(min=0.1, max=0.9), "p"),               # "p", probability of being prescribed antibiotics
    c("qunif", list(min=3, max=7), "short_dur"),           # "short_dur", mean short duration of antibiotics (normal distribution)
    c("qunif", list(min=10, max=21), "long_dur"),         # "long_dur", mean long duration of antibiotics (normal distribution)
    c("qunif", list(min=1, max=5), "sdDur")               # "sdDur", standard deviation of duration of antibiotics 
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
LHS.simple <- LHS(modelRun.simple, factors, N=3000, q, q.arg, nboot=10, cl=cl) #N is the size of the hypercube
# print elapsed time
new <- Sys.time() - old # calculate difference
print(new) # print in nice format


old <- Sys.time() # get start time
check.LHS.simple <- LHS(modelRun.simple, factors, 2000, q, q.arg, nboot=10, cl=cl)
new <- Sys.time() - old # calculate difference
print(new) # print in nice format


results.simple <- get.results(LHS.simple)
# Save run to disk
image_name <- paste0("LHS_", model, "_", format(Sys.time(), "%d%b%Y_%H%M%Z"))
save.image(paste0("./runs/", image_name, ".Rdata"))

stopCluster(cl)

##################################### Display results ########################################
load('LHS_simple_05Dec2018_0211GMT.Rdata')

#Plot findings: 
#1. empirical cumulative distribution function used to illustrate the distribution of the model results
plotecdf(LHS.simple) #outcome has a high probability in the steepest parts of the graph 

#2. scatterplot of the result as a function of each parameter: distribution of values returned by the model 
# in the parameter space sampled by the hypercube and how sensible are these model responses to the variation of each parameter.
plotscatter(LHS.simple)

#3. partial correlation coefficient measures how strong are the inear associations between the result 
# and each input parameter, after removing the linear effect of the other parameters. (CI generated by boostrapping)
plotprcc(LHS.simple)

pic(LHS.simple, nboot=40) #pic (partial inclination coefficient) is the sensitivity" of the model response in respect to each parameter
#represent the beta terms in y = alpha + beta*x regressions, after removing the linear effect of the other parameters

# 4. Cobweb
outcome.df<-as.data.frame(cbind(LHS.simple$data,results.simple)) #dummy matrix with parameter values in columns and outcome in last column
names(outcome.df)<- c(factors,'outcome') #name the columns of the dummy matrix 
for (i in 1:nrow(outcome.df)) {       #label the rows of parameter values that produced top 5% of the outcomes
    if (outcome.df$outcome[i]<quantile(outcome.df$outcome,probs = 0.95)) { 
        outcome.df$top5[i] <-0 } else {
            outcome.df$top5[i] <-1
        }
}
require(plotrix) #load MASS package
blue<-alpha("lightskyblue1", alpha=0.3)
red<-alpha("red", alpha=0.6)
colors<- c(blue, red) #choose 2 colors - 1 for parameters that produced top 5% of outcomes and one for the rest
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
#check.LHS.simple <- LHS(modelRun.simple, factors, 1900, q, q.arg, nboot=20)
(mySbma <- sbma(LHS.simple, check.LHS.simple))
# value of -1 indicates complete disagreement between the runs 
# value of 1 indicated complete agreement  (>0.7 acceptable)
#caveat: if none of the model parameters is monotonically correlated with the output, the agreement between runs may stay as low as 0.2 even for very large hypercubes.
