# load libraries 
require(pse) #load pse package for Latin Hypercube
require(sensitivity) #load sensitivity package for sensitivity analysis 
require(parallel) # load parallel processing package to use multiple cores on computer (or cluster)
library(spartan) #load spartan package for AA and eFAST

################################### Consistency testing ############################################
#resource: https://cran.r-project.org/web/packages/spartan/vignettes/sensitivity_analysis.html
#data download: http://www.kieranalden.info/index.php/spartan/

setwd('/Users/moyin/Desktop/angelsfly/indiv_abxduration/')

# SAMPLE PARAMETER SPACE 
# source functions on all cores
cl <- makeCluster(detectCores()-1)
clusterCall(cl, function() {source('model_simple.R')})

#parameter space for simple model 
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
  c("qunif", list(min=3, max=50), "n.bed"),              #"n.bed", number of beds in the ward
  c("qunif", list(min=3, max=30), "mean.max.los"),       #"mean.max.los", mean of length of stay
  c("qunif", list(min=0, max=0.9), "prob_StartBact_R"),  #"prob_StartBact_R",probability of initial carriage of resistant organisms
  c("qunif", list(min=0, max=1), "prop_S_nonR"),         #"prop_S_nonR", proportion of S in the population of S and ss
  c("qunif", list(min=0, max=1), "bif"),                 #"bif", bacterial interference factor
  c("qunif", list(min=0, max=0.05), "pi_ssr"),           # "pi_ssr" probability of being transmitted r to ss (ss—> ssr)
  c("qunif", list(min=0, max=0.2), "repop.s1"),          # "repop.s1" probability of ss repopulated to S (Palleja, Nature Biology, 2018 on gut recovery ~9 months)
  c("qunif", list(min=0, max=0.005), "mu_r"),            # "mu_r", probability of decolonisation (Haggai Bar-Yoseph, JAC, 2016, decreasing colonization rates from 76.7% (95% CI=69.3%–82.8%) at 1 month to 35.2% (95% CI=28.2%–42.9%) at 12 months of follow-up)
  c("qunif", list(min=0.1, max=0.9), "abx.clear"),       # "abx.clear", probability of S becoming ss after being on antibiotics
  c("qunif", list(min=0.1, max=0.9), "p"),               # "p", probability of being prescribed antibiotics
  c("qunif", list(min=3, max=7), "short_dur"),           # "short_dur", mean short duration of antibiotics (normal distribution)
  c("qunif", list(min=8, max=20), "long_dur"),           # "long_dur", mean long duration of antibiotics (normal distribution)
  c("qunif", list(min=1, max=2), "sdDur")                # "sdDur", standard deviation of duration of antibiotics
)  

# arrange parameters in a way LHS will be happy with
q <- unlist(lapply(parameters, function(l) l[[1]]))
q.arg <- lapply(parameters, function(l) l[2:3])
factors <- unlist(lapply(parameters, function(l) l[[4]]))

# MAKE SURE the variable listing and ORDER MATCHES the variable listing input into diff_prevalence
source(paste0("model_simple.R"))
if(!(sum(factors == parameters_simple) ==  length(parameters_simple))){
  stop("Test Error: Listing of parameters in cobweb does not match parameters accepted by diff_prevalence function.")
}

# Use the LHS function to generate a hypercube 
iterationstotry= c(1, 5, 10, 20) #iterations we are going to test 
numberofrepeatsineachiteration=20 
aa_data_simple<-list() 
for (i in 1: (max(iterationstotry)*numberofrepeatsineachiteration)){
  print(paste('Calculating', i, 'in', (max(iterationstotry)*numberofrepeatsineachiteration), 'total runs'))
  old <- Sys.time() # get start time
  LHS.simple <- LHS(modelRun.simple, factors, N=200, q, q.arg, nboot=1000, cl=cl)
  samples<-LHS.simple[['data']] # save sampled parameter combinations 
  results <- get.results(LHS.simple) #save results of the simulations
  aa_data_simple[[i]]<-cbind(samples, results) #combine sampled parameter combinations and results in one file 
  new <- Sys.time() - old # calculate difference
  print(new) # print elapsed time
} 

dirtostoreAAruns='/Users/moyin/Desktop/angelsfly/indiv_abxduration/runs/AA_simple'

#store simulation results in appropriate folders 
for (i in iterationstotry){
  data=aa_data_simple
  setwd(dirtostoreAAruns)
  numbertoputinfolderi=i*numberofrepeatsineachiteration
  dir.create(as.character(i))
  for (k in 1:numberofrepeatsineachiteration){
    directory=paste0(dirtostoreAAruns,'/',i)
    setwd(directory)
    dir.create(as.character(k))
    for (g in 1:i){
      directory=paste0(dirtostoreAAruns,'/',i, '/',k)
      setwd(directory)
      dir.create(as.character(g))
      directory=paste0(dirtostoreAAruns,'/',i,'/',k, '/',g)
      setwd(directory)
      write.csv(data[[1]],'aa_data_simple.csv',row.names=FALSE)
      data=data[-1]
    }
  }
}

stopCluster(cl)

#Running Aleatory Analysis 
FILEPATH <- dirtostoreAAruns #already in dirtostoreAAruns stated above 
# Sample sizes (number of simulation replicates in each distribution) to be analysed
SAMPLESIZES <- iterationstotry
# The simulation output measures to be analysed
MEASURES <- "results"
# Number of distributions being compared. Default: 20, as performed by Read et al
NUMSUBSETSPERSAMPLESIZE <- numberofrepeatsineachiteration
# Output file name containing the simulation responses.
RESULTFILENAME <- "aa_data_simple.csv"
# Notes the column in the CSV results file where the results start.
OUTPUTFILECOLSTART <- 14
# Last column of the output measure results
OUTPUTFILECOLEND <- 14
# The A-Test value either side of 0.5 which should be considered a 'large difference'
# between two sets of results. Use of 0.23 was taken from the Vargha-Delaney publication
LARGEDIFFINDICATOR <- 0.23
# A-Test values above 0.5 (no difference) which should be considered as small,
# medium, and large differences between two result sets. Used in the graph
# summarising all sample sizes.
SMALL <- 0.56
MEDIUM <- 0.66
LARGE <- 0.73

# A summary file is created containing the median A-Test values for each sample size.
SUMMARYFILENAME <- "AA_ATestMaxAndMedians.csv"

#calculate median values for each sample size
aa_summariseReplicateRuns(FILEPATH=FILEPATH, SAMPLESIZES=SAMPLESIZES, MEASURES=MEASURES, 
                          RESULTFILENAME=RESULTFILENAME, # Output file name containing the simulation responses
                          NUMSUBSETSPERSAMPLESIZE=NUMSUBSETSPERSAMPLESIZE,
                          OUTPUTFILECOLSTART=OUTPUTFILECOLSTART, OUTPUTFILECOLEND=OUTPUTFILECOLEND, 
                          SUMMARYFILENAME=SUMMARYFILENAME) #A summary file containing the median A-Test values for each sample size

# The results of the A-Test comparisons of the twenty subsets for each sample size
# are stored within an output file. 
ATESTRESULTSFILENAME <- "AA_ATest_Scores.csv"

# calculate A-test scores - get a csv file of the A test scores and graphs output for each sample size 
a_test_results <- aa_getATestResults(FILEPATH=FILEPATH, SAMPLESIZES=SAMPLESIZES, NUMSUBSETSPERSAMPLESIZE=NUMSUBSETSPERSAMPLESIZE, 
                                     MEASURES=MEASURES, ATESTRESULTSFILENAME=ATESTRESULTSFILENAME, LARGEDIFFINDICATOR=LARGEDIFFINDICATOR, 
                                     AA_SIM_RESULTS_FILE = SUMMARYFILENAME)

# A dataframe of summary of the max and median A test scores for each sample size 
sample_summary <- aa_sampleSizeSummary(FILEPATH, SAMPLESIZES, MEASURES, SUMMARYFILENAME, ATESTRESULTS_OBJECT = a_test_results)

# Graphs, by the sample_summary object - gives a summary graph of the sample_summary object above 
# Name of the graph which summarises the analysis results for all sample sizes.
GRAPHOUTPUTFILE <- "AA_ATestMaxes.pdf"
aa_graphSampleSizeSummary(FILEPATH, MEASURES, 300, SMALL, MEDIUM, LARGE, GRAPHOUTPUTFILE, SAMPLESUMMARY_OBJECT = sample_summary)

