thresholds = seq(0.01, 1, by = 0.01)
beta = seq(1.01, 10000, by = 0.01)

D <- matrix(NA, nrow=length(beta), ncol = length(thresholds))
for(i in 1:length(beta)){
  for(j in 1:length(thresholds)){
    D[i, j] = pbeta(thresholds[j], shape1 = 1, shape2 = beta[i])
  }
}
rownames(D) = beta
colnames(D) = format(round(thresholds, 2), nsmall = 2)

betalookup <- function(prop_R_density, r_threshold){
  
  lookingforarea = 1-prop_R_density
  lookingforthreshold = format(round(r_threshold, 2), nsmall = 2)
  
  # when a high proportion of patients are R, 
  # and threshold is high, 
  # not possible to have all patients to fall within the distribution 
  possible.area = D[, which(colnames(D)==lookingforthreshold)]
  if (lookingforarea < min(possible.area)) { # if the proportion with R is too large for the threshold 
    stop ('r_thres is too large')
  } else { 
    idx = which.min(abs(lookingforarea - possible.area))
    beta = as.numeric(row.names(D)[idx])
  }
  
  return(beta)
}

##test 
# prop_R = 0.65
# r_thres = 0.3
# a=betalookup(prop_R, r_thres)
# hist(rbeta(1000, 1, a)) #should look like exponential 
# sum(rbeta(1000, 1, a) > r_thres)/1000 # should be similar to prop_R

save(list = c('D', 'betalookup'), file = 'beta_density_lookup.Rdata')

#load D
load('beta_density_lookup.Rdata')
