# load libraries 
require(pse) #load pse package for Latin Hypercube
require(sensitivity) #load sensitivity package for sensitivity analysis 
library(spartan) #load spartan package for AA and eFAST

################################### Consistency testing ############################################
#resource: https://cran.r-project.org/web/packages/spartan/vignettes/sensitivity_analysis.html
#data download: http://www.kieranalden.info/index.php/spartan/

setwd('/Users/moyin/Documents/nBox/git_projects/indiv_abxduration/')
source('model_frequency.R')

#parameters 
parameters <- list(
    c(runif(1,min=5, max=50), "n.bed"),              #n.bed; number of beds in the ward
    c(runif(1,min=3, max=20), "max.los"),       #max.los; mean of length of stay (normal distribution)
    c(runif(1,min=0.1, max=1), "p.infect"),        #probability of being prescribed narrow spectrum antibiotic
    c(runif(1,min=10, max=1000), "cum.r.1"),        #admission day when cummulative prabability of HAI requiring abx.r is 1
    c(runif(1,min=0.1, max=1), "p.r.day1"),          #probability of being prescribed broad spectrum antibiotic on day 1 of admission 
    c(runif(1,min=2, max=27), "K"),                  # gut holding capacity, on log scale, largest R number possible is exp(300) - typical colonic bacteria 10^14 number/mL content https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4991899/
    c(runif(1,min=0.1, max=1), "total_prop"),             # mean of total starting amount of e coli on log scale
    c(runif(1,min=0,max=1), "prop_R"),              # mean of starting amount of resistant gut bacteria on log scale
    c(runif(1,min=0,max=0.002), "pi_ssr"),              # pi_ssr = daily probability of transmitting resistant E coli
    c(runif(1,min=0,max=1), "r_mean"),             # r_mean= R threshold level for tranmissibility
    c(runif(1,min=0.01,max=0.1), "r_growth"),           # r_growth = growth constant for logistic growth
    c(runif(1,min=0,max=0.9), "r_thres"),             # r_thres = threshold amount of bacteria before R can be transmitted
    c(runif(1,min=0.01,max=0.1), "s_growth"),           # r_growth = growth constant for logistic growth
    c(runif(1,min=10,max=15), "abx.s"),               # abxr_killr = amount of r killed by broad spectrum abx r
    c(runif(1,min=0,max=0.00000001), "abx.r"),               # abxr_kills = amount of s killed by broad spectrum abx r
    c(runif(1,min=3, max=7), "short_dur"),           #mean short duration of narrow spectrum antibiotics (normal distribution) 
    c(runif(1,min=14, max=21), "long_dur")           #mean long duration of narrow spectrum antibiotics (normal distribution)
)

# get factor values 
values <- as.numeric(unlist(lapply(parameters, function(l) l[[1]])))
factors <- unlist(lapply(parameters, function(l) l[[2]]))

# Use the LHS function to generate a hypercube 
iterationstotry= c(1, 50, 100, 125, 150) #iterations we are going to test 
numberofrepeatsineachiteration=20 
aa_data_freq_diff<-list() 
for (i in 1: (max(iterationstotry)*numberofrepeatsineachiteration)){
  print(paste('Calculating', i, 'in', (max(iterationstotry)*numberofrepeatsineachiteration), 'total runs'))
  old <- Sys.time() # get start time
  freq <- diff_prevalence(n.bed=values[1], max.los=values[2], p.infect=values[3], cum.r.1=values[4], 
                          p.r.day1=values[5],K=values[6], total_prop=values[7], 
                          prop_R=values[8], pi_ssr=values[9], 
                          r_mean=values[10], r_growth=values[11], r_thres = values[12],
                          s_growth=values[13], abx.s=values[14], abx.r=values[15], short_dur = values[16], long_dur = values[17])
  samples<- values #sampled parameter combinations 
  results <- freq #save results of the simulations
  aa_data_freq_diff[[i]]<-matrix(c(samples, results), byrow = TRUE, ncol = 20) #combine sampled parameter combinations and results in one file 
  colnames(aa_data_freq_diff[[i]])=c(parameters_diff_prevalence_freq, "long",'short', "RThresperbed")
  new <- Sys.time() - old # calculate difference
  print(new) # print elapsed time
} 

dirtostoreAAruns='/Users/moyin/Documents/nBox/git_projects/indiv_abxduration/runs/ATest_freq/test_scenarioB/'

#store simulation results in appropriate folders 
for (i in iterationstotry){
  data=aa_data_freq_diff
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
      write.csv(data[[1]],'aa_data_freq.csv',row.names=FALSE)
      data=data[-1]
    }
  }
}

#Running Aleatory Analysis 
FILEPATH <- dirtostoreAAruns #already in dirtostoreAAruns stated above 
# Sample sizes (number of simulation replicates in each distribution) to be analysed
SAMPLESIZES <- iterationstotry
# The simulation output measures to be analysed
MEASURES <- c("long","short","RThresperbed")
# Number of distributions being compared. Default: 20, as performed by Read et al
NUMSUBSETSPERSAMPLESIZE <- numberofrepeatsineachiteration
# Output file name containing the simulation responses.
RESULTFILENAME <- "aa_data_freq.csv"
# Notes the column in the CSV results file where the results start.
OUTPUTFILECOLSTART <- 18
# Last column of the output measure results
OUTPUTFILECOLEND <- 20
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


