#!/bin/bash
# Simple wrapper to run grid_model script based on parameters
# Usage: ./run.sh antibiotic_duration antibiotic_percent


#R CMD BATCH --vanilla "--args '$(pwd) $DUR $PER'" script_inputargs_test.R script_dur${DUR}percent${PER}_$(date +%Y-%m-%d_%H%M%S).Rout
#echo R CMD BATCH --vanilla '--args '$(pwd)'' script_dur${DUR}percent${PER}.R script_dur${DUR}percent${PER}_$(date +%Y-%m-%d_%H%M%S).Rout
#R CMD BATCH --vanilla '--args '$(pwd)'' script_dur${DUR}percent${PER}.R script_dur${DUR}percent${PER}_$(date +%Y-%m-%d_%H%M%S).Rout

echo R CMD BATCH --vanilla cobweb_simple.R indiv_simple_$(date +%Y-%m-%d_%H%M%S).Rout
R CMD BATCH --vanilla cobweb_simple.R indiv_simple_$(date +%Y-%m-%d_%H%M%S).Rout