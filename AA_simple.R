#########################################################################
#######Effect of antibiotic duration on hospitalised patients############
################Determine number of iterations per run###################
#########################################################################

setwd("/Users/moyin/Documents/nBox/git_projects/indiv_abxduration/")
rm(list=ls()) # Clean working environment

source("model_simple.R")

# load libraries 
require(pse) #load pse package for Latin Hypercube
require(sensitivity) #load sensitivity package for sensitivity analysis 
require(parallel) # load parallel processing package to use multiple cores on computer (or cluster)
require(spartan) #for AA 

################################### Consistency testing ############################################
#resource: https://cran.r-project.org/web/packages/spartan/vignettes/sensitivity_analysis.html
#data download: http://www.kieranalden.info/index.php/spartan/

# SAMPLE PARAMETER SPACE 
# source functions on all cores
cl <- makeCluster(detectCores()-1)
clusterCall(cl, function() {source('model_simple.R')})

parameters <- list(
    c(runif(1,min=5, max=50), "n.bed"),         #"n.bed", number of beds in the ward
    c(runif(1,min=3, max=20), "max.los"),       #"max.los", mean of length of stay
    c(runif(1,min=0, max=1), "prop_R"),         #"prop_R",probability of initial carriage of resistant organisms
    c(runif(1,min=0, max=1), "prop_S"),        #"prop_S", proportion of S in the population of S and ss
    c(runif(1,min=0, max=1), "bif"),            #"bif", bacterial interference factor
    c(runif(1,min=0.001, max=0.3), "pi_ssr"),      # "pi_ssr" probability of being transmitted r to ss (ss—> ssr)
    c(runif(1,min=0.005, max=0.015), "repop.s"), # "repop.s1" probability of ss repopulated to S (Palleja, Nature Biology, 2018 on gut recovery ~9 months)
    c(runif(1,min=0.002, max=0.02), "mu"),      # "mu_r", probability of decolonisation (Haggai Bar-Yoseph, JAC, 2016, decreasing colonization rates from 76.7% (95% CI=69.3%–82.8%) at 1 month to 35.2% (95% CI=28.2%–42.9%) at 12 months of follow-up)
    c(runif(1,min=0.1, max=0.5), "abx.s"),      # "abx.s", probability of S becoming ss after being on narrow spectrum antibiotics
    c(runif(1,min=0.1, max=0.5), "abx.r"),      # "abx.r", probability of R becoming ss after being on broad spectrum antibiotics
    c(runif(1,min=0.1, max=1), "p.infect"),     # "p.infect", probability of being prescribed antibiotics
    c(runif(1,min=10, max=1000), "cum.r.1"),    # admission day when cummulative prabability of HAI requiring abx.r is 1
    c(runif(1,min=0.1, max=1), "p.r.day1"),     # probability of being prescribed broad spectrum antibiotic on admission 
    c(runif(1,min=3, max=7), "short_dur"),      # "short_dur", mean short duration of antibiotics (normal distribution)
    c(runif(1,min=14, max=21), "long_dur")      # "long_dur", mean long duration of antibiotics (normal distribution)
)

# get factor values 
values <- as.numeric(unlist(lapply(parameters, function(l) l[[1]])))
factors <- unlist(lapply(parameters, function(l) l[[2]]))

# Use the LHS function to generate a hypercube 
iterationstotry= c(1, 50, 100, 125) #iterations we are going to test 
numberofrepeatsineachiteration=20 
aa_data_simple_diff<-list() 
for (i in 1: (max(iterationstotry)*numberofrepeatsineachiteration)){
  print(paste('Calculating', i, 'in', (max(iterationstotry)*numberofrepeatsineachiteration), 'total runs'))
  old <- Sys.time() # get start time
  simple <- diff_prevalence(n.bed=values[1], max.los=values[2], 
                            prop_R=values[3], prop_S_nonR=values[4], 
                            bif=values[5], pi_ssr=values[6], repop.s=values[7], mu=values[8], 
                            abx.s=values[9], abx.r=values[10], p.infect=values[11], 
                            cum.r.1=values[12], p.r.day1=values[13], short_dur=values[14], long_dur = values[15])
  samples<- values #sampled parameter combinations 
  results <- simple #save results of the simulations
  aa_data_simple_diff[[i]]<-matrix(c(samples, results), byrow = TRUE, ncol = 18) #combine sampled parameter combinations and results in one file 
  colnames(aa_data_simple_diff[[i]])=c(parameters_diff_prevalence_simple, 'long','short','results')
  new <- Sys.time() - old # calculate difference
  print(new) # print elapsed time
} 

dirtostoreAAruns="/Users/moyin/Documents/nBox/git_projects/indiv_abxduration/runs/ATest_simple/test_scenarioB/"

#store simulation results in appropriate folders 
for (i in iterationstotry){
  data=aa_data_simple_diff
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
MEASURES <- c("long",'short','results')
# Number of distributions being compared. Default: 20, as performed by Read et al
NUMSUBSETSPERSAMPLESIZE <- numberofrepeatsineachiteration
# Output file name containing the simulation responses.
RESULTFILENAME <- "aa_data_simple.csv"
# Notes the column in the CSV results file where the results start.
OUTPUTFILECOLSTART <- 16
# Last column of the output measure results
OUTPUTFILECOLEND <- 18
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


