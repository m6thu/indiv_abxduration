# =========================================================================== #
# Runs STAN models based on data in metareg_define_data.R and meta0.1.stan
# =========================================================================== #

rm(list = ls())
library(rstan)

source('define_data.R')

fit <- stan(
  file = "meta0.1JA.stan",  # Input model version here 
  data = stan_data_esbl,  # named list of data defined in define_data.R
  chains = 1,             # number of Markov chains
  warmup = 1000,          # number of warmup iterations per chain
  iter = 10000,           # total number of iterations per chain
  cores = 1,              # number of cores (use one per chain)
  refresh = 1000,         # no of runs at which progress is shown
  control = list(max_treedepth = 15, adapt_delta=0.99)
)

# Print relative risk for handwashing and masks
print(fit, pars = "OR_abx_dur")
summary(fit, pars = "OR_abx_dur")$summary
traceplot(fit, pars = "OR_abx_dur")

# Save output parameters
save(fit, file = "model0.1.Rdata")


